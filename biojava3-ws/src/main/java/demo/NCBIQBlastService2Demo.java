package demo;

import static org.biojava3.ws.alignment.qblast2.enums.BlastParameter.ENTREZ_QUERY;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.ws.alignment.qblast2.NCBIQBlastAlignmentProperties2;
import org.biojava3.ws.alignment.qblast2.NCBIQBlastOutputProperties2;
import org.biojava3.ws.alignment.qblast2.NCBIQBlastService2;
import org.biojava3.ws.alignment.qblast2.enums.BlastOutputAlignmentFormat;
import org.biojava3.ws.alignment.qblast2.enums.BlastOutputFormat2;
import org.biojava3.ws.alignment.qblast2.enums.BlastProgram;

/**
 * A simple demo showing {@link NCBIQBlastService2} usage
 * 
 * @author Gediminas Rimsa
 */
public class NCBIQBlastService2Demo {
	private static final String BLAST_OUTPUT_FILE = "blastOutput.xml";

	private static final String SEQUENCE = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS";

	public static void main(String[] args) {
		try {
			NCBIQBlastService2 service = new NCBIQBlastService2();

			// set alignment properties
			NCBIQBlastAlignmentProperties2 props = new NCBIQBlastAlignmentProperties2();
			props.setBlastDatabase("swissprot");
			props.setAlignementOption(ENTREZ_QUERY, "\"serum albumin\"[Protein name] AND mammals[Organism]");
			props.setBlastProgram(BlastProgram.BLASTP);

			// set output properties
			NCBIQBlastOutputProperties2 outputProps = new NCBIQBlastOutputProperties2();
			outputProps.setOutputFormat(BlastOutputFormat2.XML);
			outputProps.setAlignmentOutputFormat(BlastOutputAlignmentFormat.PAIRWISE);

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
			} finally {
				IOUtils.close(writer);
				IOUtils.close(reader);
				service.sendDeleteRequest(rid);
			}
		} catch (Exception e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
	}
}
