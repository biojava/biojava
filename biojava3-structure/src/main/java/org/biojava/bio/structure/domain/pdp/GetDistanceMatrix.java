package org.biojava.bio.structure.domain.pdp;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;

import org.biojava.bio.structure.StructureException;


public class GetDistanceMatrix {


	/** A set of Calpha atoms that are representing the protein
	 * 
	 * @param protein
	 */
	public  PDPDistanceMatrix getDistanceMatrix(Atom[] protein) throws StructureException{
		int[][] dist = new int[protein.length+3][protein.length+3];
		int i,j;
		double d,dt1,dt2,dt3,dt4;
		int nclose=0;
		int[] iclose = new int[protein.length*protein.length];
		int[] jclose= new int[protein.length*protein.length];
		
		if(protein.length >= PDPParameters.MAXLEN) {
			System.err.println(String.format("%d protein.len > MAXLEN %d\n",protein.length,PDPParameters.MAXLEN));
			return null;
		}
		for(i=0;i<protein.length;i++) {
			for(j=i;j<protein.length;j++) { 
				dist[i][j]=0;
				dist[j][i]=0;

				d=0;

				Atom ca1 = protein[i];
				Atom ca2 = protein[j];
				Group g1 = ca1.getGroup();
				Group g2 = ca2.getGroup();

				Atom cb1 = getCBeta(g1);
				Atom cb2 = getCBeta(g2);
				boolean hasCbeta1 = cb1 != null;
				boolean hasCbeta2 = cb2 != null;

				dt1=81;
				dt2=64;
				dt3=49;
				dt4=36;
				
				if(hasCbeta1 && hasCbeta2) {
					double distance = Calc.getDistance(cb1,cb2);
					d+= distance*distance;
				}
				else if(hasCbeta1 && ! hasCbeta2) {
					double distance = 999;
					
					distance = Calc.getDistance(cb1, ca2);
					d += distance * distance;
				}
				else if(!hasCbeta1&&hasCbeta2) {
					double distance = Calc.getDistance(ca1, cb2);
					d += distance * distance;
				}
				else if( ! hasCbeta1&&!hasCbeta2) {
					double distance = Calc.getDistance(ca1, ca2);
					d += distance * distance;
				}
				
				if(d<dt1) {
					dist[i][j]=1;
					dist[j][i]=1;
					if(d<dt2) {
						dist[i][j]=2;
						dist[j][i]=2;
						if(j-i>35) {
							iclose[nclose]=i;	
							jclose[nclose]=j;	
							nclose++;
						}
						if(d<dt3) {
							dist[i][j]=4;
							dist[j][i]=4;
							if(d<dt4) {
								dist[i][j]=6;
								dist[j][i]=6;
							}
						}
					}
				}
			}
		}
		/* secondary structure interaction */
		for(i=1;i<protein.length;i++) {
			for(j=i;j<protein.length-1;j++) { 
				/* beta-sheet */
				if(dist[i][j]>=2&&j-i>5) {
					if(dist[i-1][j-1]>=2&&dist[i+1][j+1]>=2||dist[i-1][j+1]>=2&&dist[i+1][j-1]>=2) {
						dist[i][j]+=4;
						dist[j][i]+=4;
						/*
					printf("1: %d %d %d\n",i,j,dist[i][j]);
						 */
					}
					/* alpha-helices */
					else if(i>2&&j<protein.length-2) {
						if(dist[i-3][j-3]>=1&&dist[i+3][j+3]>=1||dist[i-3][j+3]>=1&&dist[i+3][j-3]>=1) {
							dist[i][j]+=4;
							dist[j][i]+=4;
							/*
						printf("3: %d %d %d\n",i,j,dist[i][j]);
							 */
						}
						else if(i>3&&j<protein.length-3) {
							if((dist[i-3][j-3]>=1||dist[i-3][j-4]>=1||dist[i-4][j-3]>=1||dist[i-4][j-4]>=1)&&
									(dist[i+4][j+4]>=1||dist[i+4][j+3]>=1||dist[i+3][j+3]>=1||dist[i+3][j+4]>=1)
									||(dist[i-4][j+4]>=1||dist[i-4][j+3]>=1||dist[i-3][j+4]>=1||dist[i-3][j+3]>=1)&&
									(dist[i+4][j-4]>=1||dist[i+4][j-3]>=1||dist[i+3][j-4]>=1||dist[i+3][j-3]>=1)) {
								dist[i][j]+=4;
								dist[j][i]+=4;
								/*
							printf("4: %d %d %d\n",i,j,dist[i][j]);
								 */
							}
						}
					}
				}
			}
		}
		
		PDPDistanceMatrix matrix = new PDPDistanceMatrix();
		
		matrix.setNclose(nclose);
		matrix.setIclose(iclose);
		matrix.setJclose(jclose);
		matrix.setDist(dist);
		return matrix;

	}



	private Atom getCBeta(Group g1) {
		Atom cb = null;

		try {
			cb = g1.getAtom("CB");
		} catch (StructureException e){
			if ( g1 instanceof AminoAcid) {
				AminoAcid aa = (AminoAcid) g1;

				try {
					cb = Calc.createVirtualCBAtom(aa);
				} catch (StructureException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
		}
		return cb;
	}



	



	
	
}
