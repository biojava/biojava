package demo;

import java.net.URL;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceAligner;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DemoAlignProteins {

	private final static Logger logger = LoggerFactory.getLogger(DemoAlignProteins.class);

	public static void main(String[] args) {

		try {
			String uniprotID1 = "P69905";
			String uniprotID2 = "P68871";

			ProteinSequence s1 = getSequenceForId(uniprotID1);
			ProteinSequence s2 = getSequenceForId(uniprotID2);

			SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getBlosum65();

			GapPenalty penalty = new SimpleGapPenalty();

			short gop = 8;
			short extend = 1;
			penalty.setOpenPenalty(gop);
			penalty.setExtensionPenalty(extend);


			PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> smithWaterman =
					Alignments.getPairwiseAligner(s1, s2, PairwiseSequenceAlignerType.LOCAL, penalty, matrix);

			SequencePair<ProteinSequence, AminoAcidCompound> pair = smithWaterman.getPair();


			logger.info(pair.toString(60));

		} catch (Exception e){
			logger.error("Exception: ", e);
		}
	}

	private static ProteinSequence getSequenceForId(String uniProtId) throws Exception {
		URL uniprotFasta = new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uniProtId));
		ProteinSequence seq = FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(uniProtId);
		logger.info("id : {} {}{}{}", uniProtId, seq, System.getProperty("line.separator"), seq.getOriginalHeader());
		
		return seq;
	}
}