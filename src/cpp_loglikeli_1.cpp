#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

//' Remove NAs from a vector
//'
//' This function is needed for the direct translation of the R functions to C++. This function is 
//' taken from StackOverflow approximately: https://stackoverflow.com/questions/19156353/remove-na-values-efficiently
//'
//' @param x a vector (possibly with NAs)
//' 
//' @export
//' @examples
//' na_omit_c(c(NA,5,10,NA,303))
// [[Rcpp::export]]
arma::vec na_omit_c(arma::vec x){
  // make placeholder vector r
  arma::vec r(x.size());
  // make counter for number of non-NAs
  int k=0;
  for (unsigned j = 0; j < x.size(); j++){
    // if not NA, add 1 to counter and set r(j) = x(j)
    if (x(j) == x(j)){
      r(j) = x(j);
      k++;
    }
  }
  // resize r to non-NA size (removing NAs)
  r.resize(k);
  // subtract one because C++ indexes from 0 while inputs from R index from 1
  r -= 1.0;
  return r;
}

//' Transforms hyperparameters to priors
//' 
//' This is the C++ version of \code{\link{thetas_to_priors}}. See there for further documentation.
//' 
//' It is often about the same speed as the R version of this function due to little computational cost.
// [[Rcpp::export]]
List thetas_to_priors_c(const arma::vec& thetas, const int n, const double thresh = 1e-3) {
  // IG priors scale vector (prior on variances)
  arma::vec b = 5.0 * exp(thetas(0)) * (1 - exp(-exp(thetas(1)) / sqrt(arma::linspace(0, n - 1, n))));
  // IG priors shape vector (prior on variances)
  arma::vec a = 6.0 * arma::ones(n);
  // temporary vector to determine number of neighbors based on threshold
  // (500 is chosen as an arbitrarily large number for finding the threshold)
  arma::vec tempor = exp(-exp(thetas(2)) * arma::linspace(1, 500, 500));
  // m denotes how many nearest neighbors (the indices)
  arma::uvec m = arma::find(tempor > thresh);
  // if m has less than 2 elements or 500 elements, set m to 2
  if ((m.n_elem + 0.0) == 500) {
    m = m.head(2);
  } else if ((m.n_elem + 0.0) < 2) {m << 0 << 1;}
  // only consider first m elements of tempor
  arma::vec temp = tempor(m);
  // create matrix where each row is temp for the coefficient prior variances
  arma::mat g = arma::zeros(n, arma::as_scalar(temp.n_elem));
  g.each_row() += temp.t();
  // divide by mean of IG prior (variances) for simplicity of conjugate derivations
  g.each_col() /= (b / (a - 1));
  // return list of priors
  return List::create(a, b, g);
}

//' Gets posteriors of Bayesian methodology
//' 
//' This is the C++ version of \code{\link{get_posts}}. See there for further documentation.
//' 
// [[Rcpp::export]]
List get_posts_c(const arma::mat& datum, const arma::vec& a, const arma::vec& b,
                  const arma::mat& g, const arma::mat& NNarray) {
  // n2, N, m as usual: number of locations, number of replications and number of neighbors
  int n2 = arma::as_scalar(NNarray.n_rows);
  int N = arma::as_scalar(datum.n_rows);
  int m = arma::as_scalar(g.n_cols);
  // initialize posteriors of correct sizes
  // IG posterior(a_post, b_post) on variance of regressions
  arma::vec a_post = arma::zeros(n2);
  arma::vec b_post = arma::zeros(n2);
  // normal mean (muhat_post) and variances (G_post where each slice is for one regression) of regression coefficients
  arma::mat muhat_post = arma::zeros(n2, m);
  arma::cube G_post(m, m, n2, fill::zeros);
  // a and a_post are constant in our set-up
  a_post.fill(a(0) + N / 2.0);
  b_post(0) = b(0) + arma::as_scalar(datum.col(0).t() * datum.col(0)) / 2.0;
  // loop through all regressions
  for (int i = 1; i < n2; i++) {
    // get neighbor indices
    arma::vec gind = na_omit_c(NNarray.row(i).head(m).t());
    // nn: number of neighbors
    int nn = arma::as_scalar(gind.n_elem);
    // set-up regression as Yi ~ Xi
    arma::mat xi = -datum.cols(conv_to<uvec>::from(gind));
    arma::vec yi = datum.col(i);
    // get inverse of G_post in closed form
    arma::mat Ginv = xi.t() * xi + diagmat(1 / g.row(i).head(nn));
    // get Cholesky
    arma::mat Ginv_chol = arma::chol(Ginv);
    // tryCATCH version of this double solve to get the posterior mean of the coefficients
    // arma::vec muhat = solve(Ginv_chol, solve(Ginv_chol.t(), xi.t() * yi));
    arma::vec temp_mu, muhat;
    // see if solve(Ginv_chol.t(), xi.t() * yi) works
    bool status = solve(temp_mu, Ginv_chol.t(), xi.t() * yi);
    // if solve fails, use generalized inverse to calculate G*X'*y
    if (! status){
      muhat = pinv(Ginv) * (xi.t() * yi);
    } else {
      // see if solve(Ginv_chol, temp_mu) fails
      bool status2 = solve(muhat, Ginv_chol, temp_mu);
      // if this solve fails, use generalized inverse to calculate G*X'*y
      if (!status2){
        muhat = pinv(Ginv) * (xi.t() * yi);
      }
    }
    // fill in posteriors for returning
    muhat_post.row(i).head(nn) = muhat.t();
    G_post(span(0, nn - 1), span(0, nn - 1), span(i, i)) = Ginv_chol;
    // calculate b_post(i)
    b_post(i) = b(i) + arma::as_scalar(yi.t() * yi - muhat.t() * Ginv * muhat) / 2.0;
  }
  // return the posteriors
  return List::create(a_post, b_post, muhat_post, G_post);
}

// [[Rcpp::export]]
arma::sp_mat samp_posts_c(List posts, const arma::mat& NNarray){
  // n2: number of locations
  int n2 = arma::as_scalar(NNarray.n_rows);
  // get posteriors from list into vectors/matrices
  arma::vec ap = posts[0];
  arma::vec bp = posts[1];
  arma::mat mup = posts[2];
  // m: number of neighbors
  int m = mup.n_cols;
  // create sparse matrix
  arma::sp_mat uhat(n2, n2);
  // fill in first value with posterior mean of 1/sqrt(error)
  uhat(0,0) = (1.0 / sqrt(bp(0))) * exp(lgamma((2.0 * ap(0) + 1.0) / 2.0) - lgamma(ap(0)));
  for (int i = 1; i < n2; i++) {
    // get neighbor indices
    arma::vec gind = na_omit_c(NNarray.row(i).head(m).t());
    // get 1/sqrt(error) posterior mean and put it on the diagonal
    double tempd = 1.0 / sqrt(bp(i)) * exp(lgamma((2.0 * ap(i) + 1.0) / 2.0) - lgamma(ap(i)));
    uhat(i,i) = tempd;
    // scale the coefficients and fill them into the sparse matrix accordingly
    for (int j = 0; j < gind.n_rows + 0.0; j++){
      uhat(gind(j), i) = mup(i, j) * tempd;
    }
  }
  return uhat;
}

// [[Rcpp::export]]
double minus_loglikeli_c(const arma::vec& thetas, const arma::mat& datum, const arma::mat& NNarray){
  // get number of points n2 and number of repetitions per point N
  int n2 = arma::as_scalar(NNarray.n_rows);
  double N = arma::as_scalar(datum.n_rows);
  // IG posterior shape
  double ap = 6.0 + N / 2.0;
  // get priors from the function, turn it into matrices/vectors instead of List
  List temp_priors = thetas_to_priors_c(thetas, n2);
  arma::mat g = temp_priors[2];
  arma::vec b = temp_priors[1];
  // get m (number of neighbors), as either the maximum allowed during construction or from the prior
  double m = min(size(g)(1) + 0.0, size(NNarray)(1) + 0.0);
  // create ll for the integrated likelihood and calculate the first element to sum
  double ll = 6.0 * log(b(0)) - ap * log(b(0) + arma::accu(pow(datum.col(0), 2)) / 2);
  for (int i  = 1; i < n2; i++){
    // get m nearest neighbors
    arma::vec gind = na_omit_c(NNarray.row(i).head(m).t());
    // get Xi for regression
    arma::mat xi = -datum.cols(conv_to<uvec>::from(gind));
    // get yi for regression (response)
    arma::vec yi = datum.col(i);
    // get inverse of posterior G (scaled coefficient variances)
    arma::mat Ginv = xi.t() * xi + diagmat(1 / g.row(i).head(gind.n_elem - 0.0));
    arma::mat Ginv_chol = arma::chol(Ginv);
    // tryCATCH version of this double solve to get the posterior mean of the coefficients
    // arma::vec muhat = solve(Ginv_chol, solve(Ginv_chol.t(), xi.t() * yi));
    arma::vec temp_mu, muhat;
    // see if solve(Ginv_chol.t(), xi.t() * yi) works
    bool status = solve(temp_mu, Ginv_chol.t(), xi.t() * yi);
    // if solve fails, use generalized inverse to calculate G*X'*y
    if (! status){
      muhat = pinv(Ginv) * (xi.t() * yi);
    } else {
      // see if solve(Ginv_chol, temp_mu) fails
      bool status2 = solve(muhat, Ginv_chol, temp_mu);
      // if this solve fails, use generalized inverse to calculate G*X'*y
      if (!status2){
        muhat = pinv(Ginv) * (xi.t() * yi);
      }
    }
    // get IG posterior scale as b + (y'y - muhat' Ginvv muhat)/2
    double b_post = b(i) + (accu(pow(yi, 2)) - as_scalar(muhat.t() * Ginv * muhat)) / 2.0;
    // negative of calculate (log) determinant term sqrt(|G| / |g|)
    double ldet = 0.5 * (2.0 * accu(log(Ginv_chol.diag())) + accu(log(g.row(i).head(gind.n_elem - 0.0))));
    // calculate (log) last term b^a / b_post^a_post
    double lb = 6.0 * log(b(i)) - ap * log(b_post);
    // sum with previous log-likelihood
    ll +=  lb - ldet;
  }
  ll *= -1;
  return ll;
}