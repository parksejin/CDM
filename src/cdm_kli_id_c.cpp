//  Code created: 2014-01-24 17:43:40


// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP cdm_kli_id_C( SEXP pjk_, SEXP skillclasses) ;
}

// definition

SEXP cdm_kli_id_C( SEXP pjk_, SEXP skillclasses ){
BEGIN_RCPP
  
     /////////////////////////////////////  
     // INPUT  
     Rcpp::NumericMatrix pjk(pjk_);  
     Rcpp::NumericMatrix sc(skillclasses) ;   
       
     int I= pjk.nrow();  // number of items  
     int TP = sc.nrow() ; // number of skill classes  
     double K = sc.ncol() ; // number of attributes   
       
     double tmp1 =0 ;  
     double tmp2 = 0 ;  
     double sum_hdist = 0 ;  
//     double K2 = pow(2,K) ;  
       
     // create item wise matrices of KLI entries  
     Rcpp::NumericMatrix kli(TP,TP*I) ;  
       
     // global item discrimination  
     Rcpp::NumericVector glob_item(I) ;  
     // attribute-specific item discrimination  
     Rcpp::NumericMatrix attr_item(I,K) ;  
     Rcpp::NumericMatrix attr_item_count(I,K) ;  
       
     //***********************  
     // create distance function between skill classes  
     Rcpp::NumericMatrix hdist(TP,TP) ;  
     for (int tt=0;tt<TP;tt++){  
     for (int uu=tt+1;uu<TP;uu++){  
     	tmp1 = 0 ;  
     	for (int kk=0;kk<K;kk++){  
     		if ( sc(tt,kk) != sc(uu,kk) ){   
     			tmp1 ++ ;  
     				}  
     			}  
     	hdist(tt,uu) = tmp1 ;  
     	hdist(uu,tt) = tmp1 ;  
     	sum_hdist += 2 * tmp1 ;  
     		}  
     	}		  
       
     // index scripting of arrays  
     // [ii , jj , kk ] ~ [ ii , jj + J*kk ]   
     // [ I , J , K ]  
       
     /// compute entries of the Kullback Leibler information matrix	  
       
       
     for (int ii=0;ii<I;ii++){  
     // int ii = 0 ; // item ii  
       
     tmp2 = 0 ;  
       
     for (int tt=0;tt<TP;tt++){  
     for (int uu=0;uu<TP;uu++){  
     	// int tt = 0 ; // skill class tt  
     	// int uu = 5 ; // skill class uu  
     	kli( tt , uu + TP * ii ) = pjk( ii , tt ) * log( pjk( ii , tt ) / pjk( ii , uu ) ) +   
     		 pjk( ii , tt+TP ) * log( pjk( ii , tt+TP ) / pjk( ii , uu+TP ) ) ;  
             tmp2 += kli( tt , uu + TP*ii ) * hdist( tt , uu ) ;  
             for (int aa=0;aa<K;aa++){  
     		if ( ( sc(uu,aa) != sc(tt,aa) ) & ( hdist(uu,tt) == 1)  ){  
     			attr_item(ii,aa) += kli( tt , uu + TP*ii ) ;  
     			attr_item_count(ii,aa) ++  ;			  
     					}  
     				} // end aa  
     			}   // end uu  
     		}  // end tt  
       
     glob_item[ii] = tmp2 / sum_hdist ;   
     for (int aa=0;aa<K;aa++){  
     	attr_item(ii,aa) = attr_item(ii,aa) / attr_item_count(ii,aa) ;  
     		}  
     }  
       
     ///////////////////////////////////////  
     /// OUTPUT                  
       
       
     return List::create(_["pjk"]=pjk ,   
     		_["skillclasses"]= sc ,  
     		_["kli"] = kli , _["hdist"] = hdist ,   
     //		_["sum_hdist"] = sum_hdist ,   
     		_["I"]=I , _["TP"] = TP , _["K"] = K,  
     		_["glob_item"] = glob_item ,  
     		_["attr_item"] = attr_item ,  
     		_["attr_item_count"] = attr_item_count   
     			) ;  
     //	   _["matrk"]=MATRK  , _["indexmatr"]=INDEXMATR ) ;     
     // return List::create(_["yM"]=YM , _["wM"]=WM ) ;     
     
END_RCPP
}



