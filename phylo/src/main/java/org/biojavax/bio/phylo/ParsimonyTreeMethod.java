package org.biojavax.bio.phylo;
import java.io.*;
import java.lang.*;
import java.util.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.biojava.utils.process.ExternalProcess;
import org.biojavax.bio.phylo.io.nexus.*;
import org.jgrapht.*;
import org.jgrapht.generate.*;
import org.jgrapht.graph.*;


public class ParsimonyTreeMethod {
	
	
	public static void MP(TaxaBlock t, CharactersBlock ch){
	
      		int NTax = t.getDimensionsNTax();
		int NChar = ch.getDimensionsNChar();
     	 	List labels = t.getTaxLabels();
		
		String [] seq = new String[NTax];
		//WeightedGraph<String, DefaultWeightedEdge> [] jgrapht;
		
		//writing number of taxa & length of sequences to the phylip input file
		try{
			FileWriter fw = new FileWriter(new File("C:\\Program Files\\phylip3.67\\exe\\temp.txt"), true);
		  	fw.write(NTax + "     " + NChar + "\n");
			fw.close();
		}catch(IOException e){
			System.out.println("Error in Writing Temp_File(1)!");
		}


		for(int i = 0; i < NTax; i++){
			seq[i] = "";
		}

		int name_len = 0; // variable for finding the longest taxa name (for alignment)

		for (Iterator i = labels.iterator(); i.hasNext(); ) {
			
			String taxa = (String)i.next();
		            List matrix = ch.getMatrixData(taxa);
			
			if(name_len < taxa.length())
				name_len = taxa.length();
	
			for (Iterator j = matrix.iterator(); j.hasNext(); ) {                     
				Object elem = j.next();				
					
				if (elem instanceof Set) {
                           			Set data = (Set)elem;
				} else if (elem instanceof List) {
                             			List data = (List)elem;
				} else {
                              			String data = elem.toString();					  
					if(data != null && data != " ")
						seq[labels.indexOf(taxa)] += data;
				}
	 		}
		}

		//writing taxa name & sequence to the phylip input file	
		for(Iterator i = labels.iterator(); i.hasNext(); ) {

			String taxa = (String)i.next();

			try{
				FileWriter fw = new FileWriter(new File("C:\\Program Files\\phylip3.67\\exe\\temp.txt"), true);
			  	fw.write(taxa);

				for(int j = 0; j < name_len - taxa.length(); j++) 
					fw.write(" ");
				fw.write("     " + seq[labels.indexOf(taxa)] + "\n");
				fw.close();
			}catch(IOException e){
				System.out.println("Error in Writing Temp_File(2)!");
			}
		}

		try{
		  	FileWriter fw = new FileWriter(new File("C:\\Program Files\\phylip3.67\\exe\\temp.txt"), true);
			fw.write("\n \n");
			fw.close();
		}catch(IOException e){
			System.out.println("Error in Writing Temp_File(1)!");
		}


		ExternalProcess ep = new ExternalProcess();
		Object [] cmd = new Object[2];
		cmd[0] = "C:\\Program Files\\phylip3.67\\exe\\dnapars";
		cmd[1] = "Y";
		//cmd[2] = "F";
		StringWriter output = new StringWriter();
		
		//System.out.println(ep.joinCommands(cmd));
		
		try{
			try{
				try{
					try{
						try{
							ep.execute(ep.joinCommands(cmd), "C:\\Program Files\\phylip3.67\\exe\\temp.txt",  output, null);	
						} catch (IOException ie){}
					}catch (InterruptedException ite){}
				}catch(NullPointerException ne){}
			}catch(SecurityException se){}
                        }catch(IllegalArgumentException iae){}
		

	}	
}

