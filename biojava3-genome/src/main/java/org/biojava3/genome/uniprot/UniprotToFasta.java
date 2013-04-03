package org.biojava3.genome.uniprot;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.io.FastaWriterHelper;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.FastaWriterHelper;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Scooter
 */
public class UniprotToFasta {

    public static void main( String[] args ){
        try{
            String uniprotDatFileName = "uniprot_trembl_fungi.dat";
            String fastaFileName = "uniprot__trembel_fungi.faa";
            UniprotToFasta uniprotToFasta = new UniprotToFasta();
            uniprotToFasta.process(uniprotDatFileName, fastaFileName);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Convert a Uniprot sequence file to a fasta file. Allows you to download all sequence data for a species
     * and convert to fasta to be used in a blast database
     * @param uniprotDatFileName
     * @param fastaFileName
     * @throws Exception
     */

    public void process( String uniprotDatFileName,String fastaFileName ) throws Exception{

            FileReader fr = new FileReader(uniprotDatFileName);
            BufferedReader br = new BufferedReader(fr);
            String line = br.readLine();
            String id = "";
            StringBuffer sequence = new StringBuffer();
            ArrayList<ProteinSequence> seqCodingRegionsList = new ArrayList<ProteinSequence>();
            int count = 0;
            HashMap<String,String> uniqueGenes = new HashMap<String,String>();
            HashMap<String,String> uniqueSpecies = new HashMap<String,String>();
            while(line != null){
                if(line.startsWith("ID")){
                    String[] data = line.split(" ");
                    id = data[3];
                }else if(line.startsWith("SQ")){
                    line = br.readLine();
                    while(!line.startsWith("//")){

                        for(int i = 0; i < line.length(); i++){
                            char aa = line.charAt(i);
                            if((aa >= 'A' && aa <= 'Z') || (aa >= 'a' && aa <= 'z' )){
                                sequence.append(aa);
                            }
                        }
                        line = br.readLine();
                    }

                 //   System.out.println(">" + id);
                 //   System.out.println(sequence.toString());

                    ProteinSequence seq = new ProteinSequence(sequence.toString() );
                    seq.setAccession(new AccessionID(id));

                    seqCodingRegionsList.add(seq);
                    sequence = new StringBuffer();
                    count++;
                    if(count % 100 == 0)
                        System.out.println(count);
                    String[] parts = id.split("_");
                    uniqueGenes.put(parts[0], "");
                    uniqueSpecies.put(parts[1],"");
                }
                line = br.readLine();
            }
       //     System.out.println("Unique Genes=" + uniqueGenes.size());
       //     System.out.println("Unique Species=" + uniqueSpecies.size());
       //     System.out.println("Total sequences=" + seqCodingRegionsList.size());
            FastaWriterHelper.writeProteinSequence(new File(fastaFileName), seqCodingRegionsList);
            
            br.close();
            fr.close();

      //      System.out.println(uniqueGenes.keySet());
      //      System.out.println("====================");
      //      System.out.println(uniqueSpecies.keySet());


    }

}
