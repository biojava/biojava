package demo;

import static org.biojava3.ws.alignment.qblast.BlastAlignmentParameterEnum.ENTREZ_QUERY;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.ws.alignment.qblast.BlastProgramEnum;
import org.biojava3.ws.alignment.qblast.NCBIQBlastAlignmentProperties;
import org.biojava3.ws.alignment.qblast.NCBIQBlastOutputProperties;
import org.biojava3.ws.alignment.qblast.NCBIQBlastService;

/**
 * A simple demo showing {@link NCBIQBlastService} usage
 * 
 * @author Gediminas Rimsa
 */
public class NCBIQBlastServiceDemo {
	private static final String BLAST_OUTPUT_FILE = "blastOutput.xml";

	private static final String SEQUENCE = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS";

	public static void main(String[] args) {
		NCBIQBlastService service = new NCBIQBlastService();

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
