package org.biojavax.bio.phylo;

 /*
  *   MultipleHitCorrection methods for phylogeny inference
  *
  *   @author Bohyun Lee
  */
public class MultipleHitCorrection {
	 

       /*		calculate distance between two sequences (pairwise comparison) based on Jukes-Cantor model
	  *
	  *		@param taxa1
	  *				first sequence 
	  *
	  *		@param taxa2
	  *				second sequnce 
	  *
	  *		@returns	the calculated number in double type
	  */
	 public static double JukesCantor(String taxa1, String taxa2){
			
		taxa1 = taxa1.replace(" ", "");
		taxa2 = taxa2.replace(" ", "");
		
		int length = taxa1.length();
		
		if(length == taxa2.length()){
			//only if sequence lengths are the same, run the JC method 
		
			double counter = 0.0;	

			//for every single base pairs
			for( int i = 0 ; i < length; i++){
				//compare and increase the counter when it is not identical
				if(taxa1.charAt(i) != taxa2.charAt(i))
					counter++;
			}
							
			//calculate proportion of mismatch in the sequence 
			//and, it will be used as the probability of those two taxa which will have diff. base pair at any given site
			double p = counter/ (double) length;	
			
			//calculate evolutionary distance between them (by the formula) and return it
			return (-0.75 * Math.log(1.0-(4.0/3.0)*p));
		}else{
			System.out.println("Error: Sequence Length dose not match!\n");
			return 0.0;
		}
	}	
	
	 /*		calculate distance between two sequences (pairwise comparison) based on kimura's-2parameter model
	  *
	  *		@param taxa1
	  *				first sequence 
	  *
	  *		@param taxa2
	  *				second sequnce 
	  *
	  *		@returns	the calculated number in double type
	  */
	public static double KimuraTwoParameter(String taxa1, String taxa2){
		
		taxa1 = taxa1.replace(" ","");
		taxa2 = taxa2.replace(" ","");

		int length = taxa1.length();

		if(length == taxa2.length()){
		
			double counter1 = 0.0;
			double counter2 = 0.0;

			for( int i = 0; i < length; i++){
				
				//if two taxa have diff. base-pair at a site
				if(taxa1.charAt(i) != taxa2.charAt(i)){
					
					if((taxa1.charAt(i) == 'A' && taxa2.charAt(i) == 'G') || (taxa1.charAt(i) == 'G' && taxa2.charAt(i) == 'A')){
						
						//see if it is a transition between A and G, and if so increase counter1
						counter1++;
					}else if((taxa1.charAt(i) == 'T' && taxa2.charAt(i) == 'C') || (taxa1.charAt(i) == 'C' && taxa2.charAt(i) == 'T')){
						
						//see if it is a transition between C and T, and if so increase counter1
						counter1++;
					}else{

						//if it is not transition, then increase counter2 for the transversion
						counter2++;
					}
				}
			}	

			//calculate p and q, based on counter 1 & counter 2
			double p = counter1 / (double) length;
			double q = counter2 / (double) length;

			//calculate the distance (by formula) and return it.
			return ( (0.5)*Math.log(1.0/(1.0 - 2.0*p - q)) + (0.25)*Math.log(1.0/(1.0 - 2.0*q)));	
		}else{
			System.out.println("Error: Sequence Length dose not match!\n");
			return 0.0;
		}
	}

}

