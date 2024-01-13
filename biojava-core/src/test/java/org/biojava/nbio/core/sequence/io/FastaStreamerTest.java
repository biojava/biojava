package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Test the functionality of the {@link FastaStreamer} code
 */
public class FastaStreamerTest {

	@Test
	public void stream() throws IOException {
		String file = this.getClass().getResource("PF00104_small.fasta.gz").getFile();
		Path path = Paths.get(file);
		List<ProteinSequence> sequences;

		sequences = FastaStreamer.from(path).stream().collect(Collectors.toList());
		Assert.assertEquals("Count", 283, sequences.size());

		ProteinSequence sequence;
		sequence = sequences.get(0);
		Assert.assertEquals("A2D504_ATEGE/1-46", sequence.getOriginalHeader());
		sequence = sequences.get(sequences.size()-1);
		Assert.assertEquals("Q98SJ1_CHICK/15-61", sequence.getOriginalHeader());

		sequences = FastaStreamer.from(path)
				.batchSize(2) // Ensure there isn't an edge condition loading the next buffer
				.stream()
				.collect(Collectors.toList());
		Assert.assertEquals("Count", 283, sequences.size());
	}
}
