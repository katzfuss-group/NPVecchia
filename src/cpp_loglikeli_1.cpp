#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

// taken from stack overflow approximately: https://stackoverflow.com/questions/19156353/remove-na-values-efficiently
// [[Rcpp::export]]
arma::vec na_omitc(arma::vec x){
  arma::vec r(x.size());
  int k=0;
  for (unsigned j = 0; j < x.size(); j++){
    if (x(j) == x(j)){
      r(j) = x(j);
      k++;
    }
  }
  r.resize(k);
  r -= 1.0;
  return r;
}

// [[Rcpp::export]]
List thetas_to_priors_c(const arma::vec& thetas, const int n2, const double thresh = 1e-3) {
  arma::vec b = 5.0*exp(thetas(0))*(1 - exp(-exp(thetas(1)) / sqrt(arma::linspace(0,n2-1,n2))));
  arma::vec a = 6.0*arma::ones(n2);
  arma::vec tempor = exp(-exp(thetas(2))*arma::linspace(1,500,500));
  arma::uvec m = arma::find(tempor > thresh);
  if ((m.n_elem + 0.0) == 500) {
    m = m.head(2);
  } else if ((m.n_elem + 0.0) < 2) {m << 0 << 1;}
  arma::vec temp = tempor(m);
  arma::mat g = arma::zeros(n2, arma::as_scalar(temp.n_elem));
  g.each_row() += temp.t();
  g.each_col() /= (b/(a-1));
  return List::create(a, b, g);
}

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
  a_post.fill(a(0) + N/2.0);
  b_post(0) = b(0) + arma::as_scalar(datum.col(0).t() * datum.col(0)) / 2.0;
  // loop through all regressions
  for (int i = 1; i < n2; i++) {
    // get neighbor indices
    arma::vec gind = na_omitc(NNarray.row(i).head(m).t());
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
    bool status = solve(temp_mu, Ginv_chol.t(),xi.t() * yi);
    // if solve fails, use generalized inverse to calculate G*X'*y
    if (! status){
      muhat = pinv(Ginv) * (xi.t() * yi);
    } else {
      // see if solve(Ginv_chol, temp_mu) fails
      bool status2 = solve(muhat, Ginv_chol, temp_mu);
      // if this solve fails, use generalized inverse to calculate G*X'*y
      if (!status2){
        muhat = pinv(Ginv)*(xi.t() * yi);
      }
    }
    // fill in posteriors for returning
    muhat_post.row(i).head(nn) = muhat.t();
    G_post(span(0, nn-1), span(0, nn-1), span(i,i)) = Ginv_chol;
    // calculate b_post(i)
    b_post(i) = b(i) + arma::as_scalar(yi.t() * yi - muhat.t() * Ginv * muhat) / 2.0;
  }
  // return the posteriors
  return List::create(a_post, b_post, muhat_post, G_post);
}

// [[Rcpp::export]]
arma::sp_mat samp_posts_c(List posts, const arma::mat& NNarray){
	int n2 = arma::as_scalar(NNarray.n_rows);
	arma::vec ap = posts[0];
	arma::vec bp = posts[1];
	arma::mat mup = posts[2];
	int m = mup.n_cols;
	arma::sp_mat uhat(n2,n2);
	uhat(0,0) = (1.0 / sqrt(bp(0))) * exp(lgamma((2.0*ap(0)+1.0)/2.0) - lgamma(ap(0)));
	for (int i = 1; i < n2; i++) {
		arma::vec gind = na_omitc(NNarray.row(i).head(m).t());
		double tempd = 1.0 / sqrt(bp(i)) * exp(lgamma((2.0*ap(i)+1.0)/2.0) - lgamma(ap(i)));
		uhat(i,i) = tempd;
		// arma::vec x = mup.row(i).head(m);
		// x *= tempd;
		for (int j = 0; j < gind.n_rows + 0.0; j++){
		  uhat(gind(j),i) = mup(i, j)*tempd;
		}
		// uhat.rows(conv_to<uvec>::from(gind)).col(i) = x;
		//arma::umat locs(2,nn);
		//locs.row(0) = conv_to<urowvec>::from(gind);
		//locs.row(1).fill(i);
		//uhat(m,n,x);
	}
	return uhat;
}

// [[Rcpp::export]]
double minus_loglikeli_c(const arma::vec& thetas, const arma::mat& datum, const arma::mat& NNarray, const double N){
  int n2 = arma::as_scalar(NNarray.n_rows);
  double ap = 6.0 + N/2.0;
  List temp_priors = thetas_to_priors_c(thetas,n2);
  arma::mat g = temp_priors[2];
  arma::vec b = temp_priors[1];
  double m = min(size(g)(1) + 0.0, size(NNarray)(1) + 0.0);
  double ll = 6.0*log(b(0)) - ap*log(b(0) + arma::accu(pow(datum.col(0), 2))/2);
  for (int i  = 1; i < n2; i++){
    arma::vec gind = na_omitc(NNarray.row(i).head(m).t());
    arma::mat xi = -datum.cols(conv_to<uvec>::from(gind));
    arma::vec yi = datum.col(i);
    arma::mat Ginv = xi.t() * xi + diagmat(1 / g.row(i).head(gind.n_elem - 0.0));
    arma::mat Ginv_chol = arma::chol(Ginv);
    // arma::vec muhat = solve(Ginv_chol, solve(Ginv_chol.t(), xi.t() * yi));
    // tryCATCH IF GINV_CHOL IS NOT INVERTIBLE!!!!
    arma::vec temp_mu, muhat;
    bool status = solve(temp_mu, Ginv_chol.t(),xi.t() * yi);
    if (! status){
      muhat = pinv(Ginv) * (xi.t() * yi);
    } else {
      bool status2 = solve(muhat, Ginv_chol, temp_mu);
      if (!status2){
        muhat = pinv(Ginv)*(xi.t() * yi);
      }
    }
    double b_post = b(i) + (accu(pow(yi,2)) - as_scalar(muhat.t() * Ginv * muhat))/2.0;
    double ldet = 0.5 * (2.0 * accu(log(Ginv_chol.diag())) + accu(log(g.row(i).head(gind.n_elem - 0.0))));
    double lb = 6.0 * log(b(i)) - ap * log(b_post);
    ll +=  lb - ldet;
  }
  ll *= -1;
  return ll;
}

					
// [[Rcpp::export]]
double minus_loglikeli_c2(const arma::vec& thetas, const arma::mat& datum, const arma::mat& NNarray){
  int n2 = arma::as_scalar(NNarray.n_rows);
  double N = arma::as_scalar(datum.n_rows);
  double ap = 6.0 + N/2.0;
  List temp_priors = thetas_to_priors_c(thetas,n2);
  arma::mat g = temp_priors[2];
  arma::vec b = temp_priors[1];
  double m = min(size(g)(1) + 0.0, size(NNarray)(1) + 0.0);
  double ll = 6.0*log(b(0)) - ap*log(b(0) + arma::accu(pow(datum.col(0), 2))/2);
  for (int i  = 1; i < n2; i++){
    arma::vec gind = na_omitc(NNarray.row(i).head(m).t());
    arma::mat xi = -datum.cols(conv_to<uvec>::from(gind));
    arma::vec yi = datum.col(i);
    arma::mat Ginv = xi.t() * xi + diagmat(1 / g.row(i).head(gind.n_elem - 0.0));
    arma::mat Ginv_chol = arma::chol(Ginv);
    // arma::vec muhat = solve(Ginv_chol, solve(Ginv_chol.t(), xi.t() * yi));
    // tryCATCH IF GINV_CHOL IS NOT INVERTIBLE!!!!
    arma::vec temp_mu, muhat;
    bool status = solve(temp_mu, Ginv_chol.t(),xi.t() * yi);
    if (! status){
      muhat = pinv(Ginv) * (xi.t() * yi);
    } else {
      bool status2 = solve(muhat, Ginv_chol, temp_mu);
      if (!status2){
        muhat = pinv(Ginv)*(xi.t() * yi);
      }
    }
    double b_post = b(i) + (accu(pow(yi,2)) - as_scalar(muhat.t() * Ginv * muhat))/2.0;
    double ldet = 0.5 * (2.0 * accu(log(Ginv_chol.diag())) + accu(log(g.row(i).head(gind.n_elem - 0.0))));
    double lb = 6.0 * log(b(i)) - ap * log(b_post);
    ll +=  lb - ldet;
  }
  ll *= -1;
  return ll;
}