The matlab function run16QAP runs two sets of 16 test problems
from the QAPLIB test set that are used to benchmark the the
Frank-Wolf approach to Quadratic Assignment Problems (QAPs).

Here are the two variations on calling run16QAPs and
the output you should observe:

run16QAPs('path16');
     Problem          FW1          FW2         FW10
      chr12c        13072        13072        13072 
      chr15a        26492        17272        14214 
      chr15c        11936        11936        11936 
      chr20b         4060         3428         2786 
      chr22b         8420         7876         7314 
      esc16b          320          294          294 
       rou12       253420       238134       238134 
       rou15       371458       371458       371458 
       rou20       743884       743884       737310 
      tai10a       157954       148970       135828 
      tai15a       397376       397376       391522 
      tai17a       529134       511574       507646 
      tai20a       736140       721540       721540 
      tai30a      1894640      1894640      1844636 
      tai35a      2460940      2460940      2460940 
      tai40a      3225948      3194826      3194826 


run16QAPs('lipa16');
     Problem          FW1          FW2        FW10 
     lipa20a         3791         3779         3779 
     lipa20b        27076        27076        27076 
     lipa30a        13571        13474        13449 
     lipa30b       151426       151426       151426 
     lipa40a        32109        32109        32094 
     lipa40b       476581       476581       476581 
     lipa50a        62962        62962        62906 
     lipa50b      1210244      1210244      1210244 
     lipa60a       108488       108488       108488 
     lipa60b      2520135      2520135      2520135 
     lipa70a       171820       171785       171611 
     lipa70b      4603200      4603200      4603200 
     lipa80a       256073       255779       255779 
     lipa80b      7763962      7763962      7763962 
     lipa90a       363937       363937       363884 
     lipa90b     12490441     12490441     12490441 
     
     
For questions or comments contact     
John M. Conroy, IDA Center for Computing Sciences
conroyjohnm@gmail.com

