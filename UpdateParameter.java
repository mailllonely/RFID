package org;

import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.linear.*;


public class UpdateParameter {
	
 public static double[][] gibbs(int N,int thin,int seed)
    {
    DoubleRandomEngine rngEngine=new DoubleMersenneTwister(seed);
    Normal rngN=new Normal(0.0,1.0,rngEngine);
    Gamma rngG=new Gamma(1.0,1.0,rngEngine);
    double x=0,y=0;
    double[][] mat=new double[2][N];
    for (int i=0;i<N;i++) {
        for (int j=0;j<thin;j++) {
        x=rngG.nextDouble(3.0,y*y+4);
        y=rngN.nextDouble(1.0/(x+1),1.0/Math.sqrt(x+1));
        }
        mat[0][i]=x; mat[1][i]=y;
    }
    return mat;
    }

	
	void update_z(int [][]Counts,double [][]Phi,double [][]Theta,double []Alpha)
	{
		RealMatrix phi = new Array2DRowRealMatrix(Phi);
	    RealMatrix theta = new Array2DRowRealMatrix(Theta);
	    RealVector alpha=new ArrayRealVector(Alpha);
	    RealMatrix phitheta=phi.multiply(theta);
		int G=phitheta.getRowDimension();
		int S=phitheta.getColumnDimension();
		for(int i=0;i<G;i++)
			for(int j=0;j<S;j++)
			{
				{
				  
					BetaDistribution beta = new BetaDistribution(alpha.getEntry(i), phitheta.getEntry(i, j));
					double x = Math.random();
			        double  b = beta.inverseCumulativeProbability(x);
			        int c=Counts[2][3];
					BinomialDistribution binom= new BinomialDistribution(c,b);
					x = Math.random();
					int s=binom.inverseCumulativeProbability(x);
				     //           z[i,j] = rbinom(1,counts[i,j],p)
				}
			}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		
		    double x;
	        int b;
	        BetaDistribution beta = new BetaDistribution(40.0, 40.0);
	        BinomialDistribution binom= new BinomialDistribution(100,0.3);
			
	        int k=0;
	        for (int i = 0; i < 100; i++) {
	            x = Math.random();
	            b = binom.inverseCumulativeProbability(x);
	            k=k+b;
	            System.out.println(b);
	        }
	        System.out.println(k);
	        

	}

}

