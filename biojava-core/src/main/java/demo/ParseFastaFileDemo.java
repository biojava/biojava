package demo;


import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.util.InputStreamProvider;


/**
 * Created by andreas on 6/17/15.
 */
public class ParseFastaFileDemo {


    public ParseFastaFileDemo(){


        }

    /** e.g. download ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
     * and pass in path to local location of file
     *
     * @param args
     */
        public static void main(String[] args) {

            int mb = 1024*1024;

            //Getting the runtime reference from system
            Runtime runtime = Runtime.getRuntime();

            System.out.println("##### Heap utilization statistics [MB] #####");

            //Print used memory
            System.out.println("Used Memory:"
                    + (runtime.totalMemory() - runtime.freeMemory()) / mb);

            //Print free memory
            System.out.println("Free Memory:"
                    + runtime.freeMemory() / mb);

            //Print total available memory
            System.out.println("Total Memory:" + runtime.totalMemory() / mb);

            //Print Maximum available memory
            System.out.println("Max Memory:" + runtime.maxMemory() / mb);


            if ( args.length < 1) {
                System.err.println("First argument needs to be path to fasta file");
                return;
            }

            File f = new File(args[0]);

            if ( ! f.exists()) {
                System.err.println("File does not exist " + args[0]);
                return;
            }

            long timeS = System.currentTimeMillis();

            try {

                // automatically uncompress files using InputStreamProvider
                InputStreamProvider isp = new InputStreamProvider();

                InputStream inStream = isp.getInputStream(f);


                FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(
                        inStream,
                        new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>(),
                        new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));

                LinkedHashMap<String, ProteinSequence> b;

                int nrSeq = 0;

                while ((b = fastaReader.process(100)) != null) {
                    for (String key : b.keySet()) {
                        nrSeq++;
                        //System.out.println(nrSeq + " : " + key + " " + b.get(key));
                        if ( nrSeq % 100000 == 0)
                            System.out.println(nrSeq );
                    }

                }
                long timeE = System.currentTimeMillis();
                System.out.println("parsed a total of " + nrSeq + " TREMBL sequences! in " + (timeE - timeS));
            } catch (Exception ex) {
                Logger.getLogger(ParseFastaFileDemo.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
}
