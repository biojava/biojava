package org.biojava.bio.structure.align;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;

/** A class that clusters alternative alignments according to their
 * similarity.
 * 
 * @author Andreas Prlic
 * @since 1.5
 * @version %I% %G%
 */
public class ClusterAltAligs {
	
	public static final int DEFAULT_CLUSTER_CUTOFF = 95;
	
	
	public static void cluster(AlternativeAlignment[] aligs ){
		cluster(aligs, DEFAULT_CLUSTER_CUTOFF);
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void cluster(AlternativeAlignment[] aligs, int cutoff){
		
		
		List alist = Arrays.asList(aligs);
		List testAligs = new ArrayList(alist);
		
		List clusters = new ArrayList();
		List excludeList = new ArrayList();
		
		// check how similar the eqrs are...
		for ( int i=0 ; i< aligs.length;i++){
			AlternativeAlignment a = aligs[i];
			if ( excludeList.contains(a)){
				continue;
			}
			int[] idxA = a.getIdx1();
			
			Iterator iter = testAligs.iterator();
			List remainList = new ArrayList();
			List currentCluster = new ArrayList();
			
			currentCluster.add( new Integer(i));
			excludeList.add(a);
			
			int j=-1;
			while (iter.hasNext()){
				j++;
				AlternativeAlignment b = (AlternativeAlignment) iter.next();
				if ( excludeList.contains(b))
					continue;
				
				int[] idxB = b.getIdx1();
				
				// compare the eqrs..
				int samepos = 0;
				
				for ( int x = 0 ; x < idxA.length ;x++){
					int p1 =idxA[x];
					for (int y =0; y< idxB.length ; y++){						
						int p2 = idxB[y];
						if ( p1 == p2){
							samepos++;
						}
					}
				}
				float perpos = (samepos / (float)idxA.length) * 100;
				//System.out.println("aa " + i + " samepos:"+ samepos + 
				//		" l1:"+ idxA.length + " l2:" + idxB.length + " perpos:" + perpos);
				
				if ( perpos > cutoff){
					currentCluster.add(new Integer(j));
					excludeList.add(b);
				} else {
					remainList.add(b);
				}
				
			}		
			clusters.add(currentCluster);
			if ( remainList.size() == 0) {
				break;
			}
		}		
		
		// now print the clusters...
		
		Iterator iter = clusters.iterator();
		int cpos = 0;
		while (iter.hasNext()){
			cpos++;
			//System.out.println("cluster "+cpos+":");
			List cluster = (List) iter.next();
			Iterator iter2 = cluster.iterator();
			while (iter2.hasNext()){
				Integer i = (Integer) iter2.next();
				
				AlternativeAlignment alig = aligs[i.intValue()];
				alig.setCluster(cpos);
				//System.out.println( " ("+ aligs[i.intValue()]+")");
				
			}
		
		}
	}	
}
