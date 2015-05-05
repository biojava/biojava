package demo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.cemc.ParallelAlignment;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AtomCache;

/**
 * Demo for calculating multiple alignments in parallel, for a MultipleAlignment seed.
 * A method converts then the alignments to the closest structure 
 * 
 * @author Aleix Lafita
 *
 */
public class ParallelAlignmentDemo {

	public static void main(String[] args) throws IOException, StructureException, InterruptedException, ExecutionException {
		
		
		
		//Set the list of structure names to align
		List<String> names = Arrays.asList("1mbc", "1hlb", "1thb.a", "1ith.a");		//globins
		AtomCache cache = new AtomCache();
		int size = names.size();
		
		//Build the atomArrays of the structures
		List<Atom[]> atomArrays = new ArrayList<Atom[]>();
		for (String name:names) atomArrays.add(cache.getAtoms(name));
		
		//List to store the all-to-all alignments. Contains the n^2 pairwise alignments as a 2D matrix indicies.
		List<List<AFPChain>> afpAlignments = new ArrayList<List<AFPChain>>();
		for (int i=0; i<size; i++){
			afpAlignments.add(new ArrayList<AFPChain>());
			for (int j=0; j<size; j++)
				afpAlignments.get(i).add(null);
		}
		
	    //Initialize the executor with a variable number of threads and the Future Object List
		ExecutorService executor = Executors.newCachedThreadPool();
	    List<Future<AFPChain>> afpFuture = new ArrayList<Future<AFPChain>>();
	    
	    //Create all the possible protein pairwise combinations and call the alignment
	  	for (int i=0; i<size; i++){	  		
	  		for (int j=i+1; j<size; j++){
	    
	  			Callable<AFPChain> worker = new ParallelAlignment(atomArrays.get(i), atomArrays.get(j), new CeMain());
	  			Future<AFPChain> submit = executor.submit(worker);
	  			afpFuture.add(submit);
	  		}
	  	} 
	  	
	  	int index = 0;  //the alignment index
	  	for (int i=0; i<size; i++){
	  		for (int j=i; j<size; j++){
	  			if (i!=j){
		  			afpAlignments.get(i).add(j,afpFuture.get(index).get());
		  			afpAlignments.get(j).add(i,afpFuture.get(index).get());
		  			index++;
	  			}
	  		}
	  	}
	    executor.shutdown();
	    
	    //Display all the alignments in the console
	    for (int i=0; i<size; i++){
	    	for (int j=0; j<size; j++){
	    		if (i>j) System.out.println(AfpChainWriter.toFatCat(afpAlignments.get(i).get(j), atomArrays.get(j), atomArrays.get(i)));
	    		else if (i<j) System.out.println(AfpChainWriter.toFatCat(afpAlignments.get(i).get(j), atomArrays.get(i), atomArrays.get(j)));
	    	}
	    }
	}
}
