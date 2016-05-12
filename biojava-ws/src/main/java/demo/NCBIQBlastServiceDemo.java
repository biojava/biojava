/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package demo;

import org.biojava.nbio.core.sequence.io.util.IOUtils;
import org.biojava.nbio.ws.alignment.qblast.BlastProgramEnum;
import org.biojava.nbio.ws.alignment.qblast.NCBIQBlastAlignmentProperties;
import org.biojava.nbio.ws.alignment.qblast.NCBIQBlastOutputProperties;
import org.biojava.nbio.ws.alignment.qblast.NCBIQBlastService;

import java.io.*;

import static org.biojava.nbio.ws.alignment.qblast.BlastAlignmentParameterEnum.ENTREZ_QUERY;

/**
 * A simple demo showing {@link NCBIQBlastService} usage
 *
 * @author Gediminas Rimsa
 */
public class NCBIQBlastServiceDemo {
	private static final String BLAST_OUTPUT_FILE = "blastOutput.xml";

	private static final String SEQUENCE = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS";

	public static void main(String[] args) {
		NCBIQBlastService service = null;
        if (args.length == 1) {
            service = new NCBIQBlastService(args[0]);
        } else {
            service = new NCBIQBlastService();
        }

		// set alignment options
		NCBIQBlastAlignmentProperties props = new NCBIQBlastAlignmentProperties();
		props.setBlastProgram(BlastProgramEnum.blastp);
		props.setBlastDatabase("swissprot");
		props.setAlignmentOption(ENTREZ_QUERY, "\"serum albumin\"[Protein name] AND mammals[Organism]");

		// set output options
		NCBIQBlastOutputProperties outputProps = new NCBIQBlastOutputProperties();

		// Example of two possible ways of setting output options (in this case, it was already set by constructor)
//		outputProps.setAlignmentNumber(100);
//		outputProps.setOutputOption(BlastOutputParameterEnum.ALIGNMENTS, "100");

		String rid = null;
		FileWriter writer = null;
		BufferedReader reader = null;
		try {
			// send blast request and save request id
			rid = service.sendAlignmentRequest(SEQUENCE, props);

			while (!service.isReady(rid)) {
				System.out.println("Waiting for results. Sleeping for 5 seconds");
				Thread.sleep(5000);
			}

			// read results when they are ready
			InputStream in = service.getAlignmentResults(rid, outputProps);
			reader = new BufferedReader(new InputStreamReader(in));

			File f = new File(BLAST_OUTPUT_FILE);
			System.out.println("Saving query results in file " + f.getAbsolutePath());
			writer = new FileWriter(f);

			String line;
			while ((line = reader.readLine()) != null) {
				writer.write(line + System.getProperty("line.separator"));
			}
		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
		} finally {
			IOUtils.close(writer);
			IOUtils.close(reader);
			service.sendDeleteRequest(rid);
		}
	}
}
