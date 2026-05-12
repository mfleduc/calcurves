###### Hypothesis testing functions
#'@title Two-sided Bayes factor test for calibration
#'@description Calculates the Bayes factor \eqn{BF_{10}} for the alternative and null models 
#'for calibration testing. The null model is that the vector \eqn{\vec{p_0}} describes the 
#'true gap probabilities, and the alternative is that it does not. Letting X be the data, n be the total observations,
#'and m_1 and m_0 be the alternative and null models, the Bayes factor is
#'\eqn{
#'BF_{10} = \frac{m_1(X)}{m_0(X)}=\frac{\Gamma(\sum_i\gamma_i)}{\Gamma(n+\sum_i\gamma_i)}\frac{\prod_i \Gamma(n_i+\gamma_i)/\Gamma(\gamma_i)}{(p_{0,i})^n_i}
#'}
#'where \eqn{\vec{\gamma}} is the vector of parameters for the Dirichlet prior distribution. This is an analogue of the frequentist likelihood ratio test, and has many interpretations. 
#'See Kass and Raftery 1995 for one potential interpretation.
#'@param p0 vector of null category probabilities. 
#'@param data vector of the same length as p0 containing the number of observed members in each category. 
#'@param dir_prior prior parameters for the dirichlet distribution. Default is all \eqn{\gamma_i=1/2} which is the Jeffreys prior. 
#'@returns a list:
#' $logfactor The natural log of the two-sided Bayes factor BF_{10}
#' $interpretation The interpretation of the Bayes factor as described in Kass and Raftery 1995
two_sided_bayes_factor = function(p0,data, dir_prior = array(0.5, dim=c(1, length(p0)))){
  stopifnot(length(p0)==length(dir_prior))
  stopifnot(length(p0)==length(data))
  ## Now do the factor checking
  m0 = sum(log(p0)*data)
  m1 = lgamma(sum(dir_prior))-lgamma(sum(data)+sum(dir_prior))+sum(lgamma(data+dir_prior)-lgamma(dir_prior))
  BF10 = (m1-m0)
  if(BF10<log(1)){
    interpretation="Supports null"
  }else if(BF10<=log(3.2)){interpretation="Barely worth mentioning"}
  else if(BF10<=log(10)){interpretation="Substantial support for H1"}
  else if(BF10<=log(100)){interpretation="Strong support for H1"}
  else{interpretation="Decisive support for H1"}
  return(list(logfactor=BF10, interpretation=interpretation))
}
#'@title Two-sided likelihood ratio test for calibration
#'@description A two-sided likelihood ratio test for model calibration. Similar in theory to the Bayes factor test, this 
#'tests the hypothesis that \eqn{H_0: p=p_0} against \eqn{H_1:p\ne p_0}. We present the G^2/q' statistic presented in Smith et al 1981 given by
#'
#'\eqn{
#'G^2=\frac{2}{q'}\sum_{i=1}^{d+1}n_i\log(n_i/(np_{0,i}))
#'}
#'where 
#'
#'\eqn{q'=1+\frac{1}{6dn}\left(\sum_i \frac{1}{p_{0,i}}-1+\frac{1}{n}\sum_i\left(\frac{1}{p_{0,i}}-\frac{1}{p_{0,i}^2}\right)\right)}
#'
#'which has an asymptotic \eqn{\chi^2} distribution with d degrees of freedom in the case where there are d+1 categories.
#'@param p0 vector, the category probabilities under the null hypothesis
#'@param data vector of observed counts in each category.
#'@param alpha the desired significance level
#'@return a list
#' $pval The calculated p-value
#' $result A boolean, TRUE if the nul is rejected, FALSE if it is not.
#'@export
adj_g2_test =  function(p0,data,alpha){
  stopifnot(length(p0)==length(data))
  d = length(p0)-1
  n=sum(data)
  gsq = 2*sum(data*log(data/p0/n))
  qprime = 1+(1/(6*d*n))*(sum(1/p0)-1+1/n*sum(1/p0-1/p0^2))
  gqsp = gsq/qprime
  critical_val = qchisq(1-alpha,d)
  result = gqsp>critical_val
  pval = 1-pchisq(gqsp,d)
  return(list(pval=pval,result=result))
}
#'@title Truncated-T test for the calibration
#'@description A two-sided calibration test for model calibration based on the test derived in Balakrishnan and Wasserman 2019.
#'If there are d+1 categories, which corresponds to d nominal coverage levels, the test statistic \eqn{T_{trunc}} is given by
#'\eqn{
#'T_{trunc}=\sum_i\frac{(n_i-np_{0,i})^2-n_i}{t_i}
#'}
#'where
#'\eqn{t_i=max(1/(d+1), p_{0,i})}.
#'The null hypothesis is rejected at the alpha level if
#'\eqn{
#'T_{trunc}>n\sqrt{\frac{2}{\alpha}\sum_i\left(\frac{p_{0,i}}{t_i}\right)^2}
#'}
#'@param p0 vector, the category probabilities under the null hypothesis
#'@param data vector of observed counts in each category.
#'@param alpha the desired significance level
#'@return a list
#' $cv The critical value 
#' $ttrunc The value of \eqn{T_{trunc}}
#' $result A boolean, TRUE if the nul is rejected, FALSE if it is not.
#'@export
t_trunc_test =  function(p0,data,alpha){
  stopifnot(length(p0)==length(data))
  d = length(p0)-1
  n=sum(data)
  ti = pmax(1/(d+1), p0)
  ttrunc = sum(((data-n*p0)^2-data)/ti)
  crit_val = n*sqrt(2/alpha)*sqrt(sum(p0^2/ti^2))
  return(list(cv=crit_val,ttrunc=ttrunc,result=(ttrunc>crit_val)))
}
