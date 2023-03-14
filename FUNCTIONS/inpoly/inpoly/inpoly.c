/*
 
 Returns points laying inside/on bondary of a given polygon.
 Mex file of inpoly function of Darren Engwirda. 
 
  Usage
  -----

  [in , bnd] = inpoly(X , node , [cnect])
 
 Inputs
 ------
 
   X                  Point to test if inside/bounds (2 x N).
   node               Node of Polygon (2 x Nnode) where node(: , i) is linked with node(: , i+1).
   cnect              ndex of connected noded (default node  = [(1:Nnode-1) , Nnode ; (2:Nnode) , 1])
 
 
 Outputs
 -------
 
   in                 Logical array of X points inside the polygon (1 x N).
   bnd                Logical array of X points laying on the boundary of the polygon (1 x N).
 
 
 To Compile
 ----------
   mex  inpoly.c
 
 Author
 ------
 
 Sébastien PARIS (sebastien.paris@lsis.org)
 

X          = [-81.8035 -81.8144 -81.8144 -81.8134;24.55 24.5509 24.5599 24.55];
node       = [-81.8162 -81.8318 -81.8096 -81.8061 -81.8162;24.5695 24.551 24.5453 24.5709 24.5695];
[in , bnd] = inpoly(X , node); 


 
Example 1
---------
 
close all, clc
 
N          = 200;
dtheta     = pi/150;
theta      = (-pi:dtheta:(pi-dtheta));
node       = [cos(theta) ; sin(theta)];
X          = 3*(rand(2 , N)-0.5);
[in , bnd] = inpoly(X , node);
 
plot(X(1 , in), X(2 , in),'b.',X(1 , ~in), X(2 , ~in),'r.' , node(1 , :) , node(2 , :))
title('Inside points (blue) & outside points (red)')
 
 
Example 2
---------
 
 
clc, close all
N          = 5000;
X          = [20*(rand(1 , N)-0.5) ; 10*(rand(1 , N))];
p1         = [
-8.9154147   1.6615920
-8.8847064   1.5538653
-8.7270571   1.4396306
-8.6780329   1.3310253
-8.2763869   1.3642085
-7.7589437   1.5521208
-7.3399811   1.5866412
-7.2064877   1.7400111
-6.9871610   1.7838988
-6.9190874   1.6751018
-6.6594362   1.8237880
-6.4705078   2.0285109
-6.3982459   2.0257761
-6.1776033   2.1236043
-5.9683636   1.9041636
-6.0500631   1.6420917
-6.2751037   1.4382319
-6.2083449   1.2768143
-6.3593873   1.1233594
-6.4379782   0.9672913
-6.4389687   0.9408290
-6.3308569   0.9103191
-5.9216680   1.1606956
-5.8853530   1.1594465
-5.5306989   0.8827500
-5.4571118   0.9068911
-5.2776995   0.8218321
-4.6882658   1.0696109
-4.4325551   1.1157840
-3.9867787   1.5284986
-3.6232783   1.5733560
-3.5156848   1.5181578
-3.1151141   1.6162773
-2.6031886   1.9253444
-2.4201615   2.0814581
-2.0576784   2.1825510
-1.9469493   2.3929843
-1.6564457   2.6016249
-1.4385855   2.8113568
-1.1857180   2.9681802
-0.8965760   3.2839699
-0.4656316   3.4939341
-0.0715676   3.6520868
 0.3221563   3.5994536
 0.6089023   3.4943945
 0.6092879   3.3885375
 0.3942456   3.3878929
 0.2866788   3.4141384
 0.0717038   3.3345153
 0.0896864   3.2286668
 0.1076576   3.1757488
-0.1436793   3.0169900
-0.3235837   2.8584538
-0.4676923   2.7529346
-0.9391942   2.0668111
-1.0495193   1.7498742
-1.0882363   1.3531376
-1.0158452   1.3262275
-0.9424102   1.4845979
-0.6152380   1.7477533
-0.5061946   1.9061816
-0.2890742   2.0115316
-0.5066670   1.7473959
-0.6160022   1.5360392
-0.3979720   1.8000363
-0.1808403   1.8525961
 0.2894343   1.7998173
 0.5431103   1.6681146
 0.9058832   1.5373304
 0.9432856   1.3258124
 1.0534256   1.1147328
 1.1641909   0.8507940
 1.3109203   0.6930545
 1.3117251   0.5871977
 1.4214690   0.5351345
 1.9317215   0.5400872
 2.2226054   0.5965222
 2.3340522   0.4391798
 2.4073552   0.4137156
 2.4449497   0.3348350
 2.5909063   0.3369578
 2.6249835   0.4962917
 2.7343488   0.4979874
 2.8454551   0.3938967
 2.9366470   0.3954205
 2.9201944   0.2892562
 2.9931921   0.2905022
 3.3135431   0.7199108
 3.3863589   0.7213257
 3.9960329   1.1315466
 4.0705567   1.0538628
 4.6490632   1.1744760
 4.9039980   1.1550887
 5.4495532   1.1450567
 5.9850466   1.4278414
 6.1653805   1.4607402
 6.3476152   1.4409063
 6.6355959   1.5048353
 6.5309084   1.3948027
 6.4623697   1.2862027
 6.4803085   0.8098819
 6.4094979   0.7542354
 6.4100090   0.7012564
 6.4882671   0.5981826
 6.7068647   0.6065329
 6.8546905   0.5593290
 6.9296957   0.5093111
 7.3625916   0.6332439
 7.4195041   0.5826584
 7.4810256   0.4262426
 7.6622000   0.4606821
 7.7339418   0.4903887
 7.6927729   0.5946095
 7.6199031   0.5913803
 7.5788199   0.6956243
 7.4649955   0.7967085
 7.4604098   0.9025555
 7.6668095   1.1767640
 7.7709460   1.2875089
 7.7298813   1.3917132
 7.6707366   1.4951115
 7.6356734   1.4670289
 7.6246141   1.3074531
 7.5532118   1.2777792
 7.5280861   1.4357493
 7.5345391   1.7011603
 7.7131429   1.7621690
 7.8603274   1.7158544
 7.9228745   1.9308907
 7.7348899   2.0813402
 7.3713821   2.1182087
 7.2698318   1.9547751
 7.1885795   2.1633926
 7.1152657   2.1868122
 7.0441691   2.1573400
 6.9336956   2.2058041
 6.9676141   2.2602054
 7.0648074   2.5292990
 7.0897106   2.7954185
 7.2917362   3.1221593
 7.3230152   3.2295634
 7.2074219   3.4101256
 7.0987913   3.4319941
 7.0562590   3.5892427
 6.8349672   3.7390777
 6.5426013   3.8863685
 6.3538119   4.1440578
 6.6493349   4.7917061
 6.5385166   4.8932842
 6.7052993   5.1650313
 6.7009764   5.2708721
 5.7441174   5.2345817
 5.5686945   5.1755508
 4.7886070   5.1509481
 4.4636817   5.3536511
 4.2075761   5.6646354
 3.6075780   7.1858876
 3.3873644   7.6576363
 3.2804384   7.7613007
 3.1762847   7.7327191
 3.0381962   7.6506117
 2.8286881   7.6467634
 2.8598071   7.8590921
 2.6845531   7.9089847
 2.5468001   7.8008409
 2.4437671   7.6933492
 2.3390419   7.6917801
 2.1983122   7.7691861
 2.1616242   7.9010276
 2.0918985   7.9000964
 2.0235167   7.7933405
 1.9540682   7.7660064
 1.5705080   7.7353016
 1.5002116   7.7875578
 1.2547448   7.9442315
 0.7309650   8.1526824
 0.4522013   8.2575018
 0.1738657   8.3098811
-0.0694533   8.5215153
-0.1737495   8.4157382
-0.2084645   8.4422444
-0.3476150   8.3630957
-0.4518238   8.3898231
-0.5558115   8.4695476
-0.5215963   8.3106439
-0.5221188   8.1518584
-0.4879606   7.9400337
-0.2092654   7.8335660
-0.2097517   7.4630660
-0.5257651   7.0403597
-0.5604465   7.1463352
-0.7716285   6.9354925
-0.9473082   6.8835005
-0.9833619   6.7249256
-1.1593398   6.6731632
-1.2287946   6.7795396
-1.1574372   6.9378046
-0.6994015   7.4115277
-0.5933112   7.7286747
-0.6973173   7.8878837
-0.9759166   7.9422776
-1.1175557   7.6256239
-1.1194069   7.3609825
-1.2951661   7.2564438
-1.5483098   6.4117610
-1.9366363   6.3100700
-1.8637764   6.5209552
-1.7582839   6.5197737
-1.7202965   6.7840348
-1.5427314   6.9939696
-1.5417155   7.0998256
-1.5762336   7.1530937
-1.9985176   6.9989869
-2.3506579   6.8978706
-2.5645144   6.6893724
-2.5351899   6.3183276
-2.7562921   5.7395724
-3.7980607   5.0191439
-4.3373509   4.7676755
-4.9834948   4.5740477
-5.9256532   4.1282314
-7.0382285   3.1643631
-7.4517045   2.7580044
-8.2033903   2.1245258
-8.7235361   1.8639192
]';
 
 
p2           =[
-0.1745905   7.6482740
-0.1048759   7.4629629
-0.0174793   7.4629295
-0.0349355   7.5687895
 0.0349007   7.7275752
 0.2267794   7.7806613
 0.3839078   7.7280333
 0.5583180   7.7550129
 0.6632222   7.7024853
 0.8201638   7.7296796
 0.7325480   7.8086476
 0.6273773   7.9405222
 0.4880418   7.9135694
 0.1393711   7.9922753
 0.0348428   7.9922181
-0.1743299   7.8864526
]';
 
 
n1         = size(p1 , 2);
c1         = [(1:n1-1) , n1 ; (2:n1) , 1];
n2         = size(p2 , 2);
c2         = [(1:n2-1) , n2 ; (2:n2) , 1];
 
 
node       = [p1 , p2];
cnect      = [c1 , c2 + n1];
 
[in , bnd] = inpoly(X , node , cnect );
 
plot(X(1 , in), X(2 , in),'b.',X(1 , ~in), X(2 , ~in),'r.' , node(1 , :) , node(2 , :))
title('Inside points (blue) & outside points (red)')
 
 */

#include <math.h>
#include "mex.h"

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

/*-------------------------------------------------------------------------------------------------------------- */
/* Function prototypes */

void qsindex( double * , int * , int , int  );

/*-------------------------------------------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray** plhs ,  int nrhs, const mxArray** prhs)
{
    double *X , *node, *cnect;
	double *Xtemp , *nodetemp;
    bool *in , *bnd;
    bool *cn , *on;
    double *x, *y;
    int *index;
    double tol, tol1 , eps = 2.2204e-016, temp , minix=10e20, maxix=-10e20, miniy=10e20, maxiy=-10e20, norminf=-10e20;
    double x1, x2 , y1, y2, xmin, xmax , ymin, ymax, XX, YY , x2x1 , y2y1;
    double lim , ub;
    int i , j, l , i2 , ind , k ;
    int n1 , n2 , lower , upper, start;
    int N, N1, Nnode, Ncnect;

    if (nrhs < 2)
	{
		mexPrintf(
			"\n"
			"\n"
			"Returns points laying inside/on bondary of a given polygon.\n"
			"Mex file of inpoly function of Darren Engwirda.\n" 
			"\n"
			"\n" 
			"Usage\n"
			"-----\n"
			"\n"
			"\n"
			"[in , bnd]        = inpoly(X , node , [cnect]);\n"
			"\n"
			"\n"
			"Inputs\n"
			"------\n"
			"\n"
			"X                  Point to test if inside/bounds (2 x N).\n"
			"node               Node of Polygon (2 x Nnode) where node(: , i) is linked with node(: , i+1).\n"
			"cnect              ndex of connected noded (default node  = [(1:Nnode-1) , Nnode ; (2:Nnode) , 1]).\n"
			"\n"
			"\n"
			"Outputs\n"
			"-------\n"
			"\n"
			"in                 Logical array of X points inside the polygon (1 x N).\n"
			"bnd                Logical array of X points laying on the boundary of the polygon (1 x N).\n"
			"\n"
			"\n"
			);
		return;


	}
      
    X         = mxGetPr(prhs[0]);  /* (2 x N) vector */
    if(mxGetM(prhs[0]) != 2)
    {
//        mxErrMsgTxt("X must be (2 x N)");   
    }
    
    N         = mxGetN(prhs[0]);
	N1        = N - 1;

	Xtemp     = (double *) malloc(2*N*sizeof(double));
	for(i = 0 ; i < 2*N ; i++)
	{
		Xtemp[i] = X[i];
	}
    
    node      = mxGetPr(prhs[1]);  /* (2 x Nnode) vector */
    if(mxGetM(prhs[1]) != 2)
    {
 //       mxErrMsgTxt("node must be (2 x Nnode)");   
    }
    
    
    Nnode     = mxGetN(prhs[1]);
	nodetemp  = (double *) malloc(2*Nnode*sizeof(double));
	for(i = 0 ; i < 2*Nnode ; i++)
	{
		nodetemp[i] = node[i];
	}
    
    if (nrhs < 3)
    {
        Ncnect   = Nnode;
        cnect    = malloc(2*Ncnect*sizeof(double));
        for (i = 0 ; i < Ncnect - 1 ; i++)
        {
            i2            = 2*i;
            cnect[i2]     = i + 1;
            cnect[i2 + 1] = i + 2;
        }
        cnect[Nnode*2-2] = Nnode;
        cnect[Nnode*2-1] = 1.0;
    }
    else
    {
        cnect  = mxGetPr(prhs[2]); /* (2 x Ncnect) vector */
        Ncnect = mxGetN(prhs[2]);
    }
    
    /*------------- Outputs -----------*/
    
    plhs[0]   = mxCreateLogicalMatrix(1 , N);
    in        = (bool *)mxGetPr(plhs[0]);

	plhs[1]   = mxCreateLogicalMatrix(1 , N);
    bnd       = (bool *)mxGetPr(plhs[1]);
      
    /* Temporal vectors */
    
    cn        = (bool *)malloc(N*sizeof(bool));
    on        = (bool *)malloc(N*sizeof(bool));
    x         = (double *)malloc(N*sizeof(double));
    y         = (double *)malloc(N*sizeof(double));
    index     = (int *)malloc(N*sizeof(int));
    
    for (i = 0; i < 2*N ; i=i+2)
    {        
        temp  = Xtemp[i];
        if(temp < minix)
        {
            minix = temp;
        }
        if(temp > maxix)
        {
            maxix = temp;
        }
        temp  = X[i + 1];
        if(temp < miniy)
        {
            miniy = temp;
        }
        if(temp > maxiy)
        {
            maxiy = temp;
        }
    }
    for (i = 0 ; i < 2*Nnode ; i = i + 2)
    {
        temp = fabs(node[i]);
        if(temp > norminf)
        {
            norminf = temp;
        }
        temp = fabs(node[i + 1]);
        if(temp > norminf)
        {
            norminf = temp;
        }
    }
    tol  = norminf*eps;
    lim  = norminf + tol;
	tol1 = tol + 1.0;
    
    if ((maxix - minix) > (maxiy - miniy))
    {
        for (i = 0; i < 2*N ; i=i+2)
        {
            temp             = Xtemp[i];
            Xtemp[i]         = Xtemp[i + 1];
            Xtemp[i + 1]     = temp;
        }
        for (i = 0 ; i < 2*Nnode ; i=i+2)
        {
            temp             = nodetemp[i];
            nodetemp[i]      = nodetemp[i + 1];
            nodetemp[i + 1]  = temp;
        }
    }
    for (i = 0 ; i< N ; i++)
    {
        cn[i]    = 0;
        on[i]    = 0;
        y[i]     = Xtemp[2*i + 1];
        index[i] = i;
    }
    
    qsindex( y , index , 0 , N1 );
    for (i = 0 ; i < N ; i++)
    {
        x[i]     = Xtemp[2*index[i]];
    }
    
    for (k = 0 ; k < Ncnect ; k++)
    {
        n1   = 2*(((int) cnect[2*k]) - 1);
        n2   = 2*(((int) cnect[2*k + 1]) - 1);
        x1   = nodetemp[n1];
        y1   = nodetemp[n1 + 1];
        x2   = nodetemp[n2];
        y2   = nodetemp[n2 + 1];
		x2x1 = x2 - x1;
		y2y1 = y2 - y1;

        if (x1 > x2)
        {
            xmin = x2;
            xmax = x1;
        }
        else
        {
            xmin = x1;
            xmax = x2;
        }
        if (y1 > y2)
        {
            ymin = y2;
            ymax = y1;
        }
        else
        {
            ymin = y1;
            ymax = y2;
        }
        if (y[0] == ymin)
        {
            start = 0;
        }
        else if (y[N1] <= ymin)
        {
            start = N1;
        }
        else
        {
            lower = 0;
            upper = N1;
            start = ((lower+upper)/2);
            for (l = 0 ; l < N ; l++)
            {
                if (y[start] < ymin)
                {
                    lower = start;
                    start = ((lower+upper)/2);
                }
                else if (y[start-1]<ymin)
                {
                    break;
                }
                else
                {
                    upper = start;
                    start = ((lower+upper)/2);
                }
            }
            start--;
        }

		start = max(0 , start);
        for (j = start ; j < N ; j++)
        {
            YY = y[j];
            if (YY <= ymax)
            {
                if (YY >= ymin)
                {
                    XX = x[j];
                    
                    if (XX >= xmin)
                    {
                        if (XX <= xmax)
                        {
                            on[j] = on[j] || ( fabs( (y2 - YY)*(x1 - XX) - (y1 - YY)*(x2 - XX) ) < tol );
                            if (YY < ymax)
                            {
                                ub = ( x2x1*(y1 - YY) - y2y1*(x1 - XX) )/( (XX - lim)*y2y1 );
                                if ( (ub > -tol) && (ub < tol1 ) )
                                {
                                    cn[j] = !cn[j];
                                }
                            }
                        }
                    }
                    else if (YY < ymax)   
                    {
                        cn[j] = !cn[j];   
                    }    
                }
            }
            else   
            {
                break;
            }   
        }    
    }
    for(i = 0 ; i < N ; i++)
    {
        ind      = index[i];
        in[ind]  = (cn[i] || on[i]);
        bnd[ind] = on[i];   
    }
    free(cn);
    free(on);
    free(x);
    free(y);
    free(index);
	free(Xtemp);
	free(nodetemp);
	if (nrhs < 3)
	{
		free(cnect);
	}
}

/* ----------------------------------------------------------------------------- */
void qsindex (double  *a, int *index , int lo, int hi)
{
/*  
  lo is the lower index, hi is the upper index
  of the region of array a that is to be sorted
*/
    int i=lo, j=hi , ind;
    double x=a[(lo+hi)/2] , h;

    /*  partition */
    do
    {    
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            h        = a[i]; 
			a[i]     = a[j]; 
			a[j]     = h;
			ind      = index[i];
			index[i] = index[j];
			index[j] = ind;
            i++; 
			j--;
        }
    }
	while (i<=j);

    /*  recursion */
    if (lo<j) qsindex(a , index , lo , j);
    if (i<hi) qsindex(a , index , i , hi);
}
/* ----------------------------------------------------------------------------- */
