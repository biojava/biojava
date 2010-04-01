/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.benchmark.metrics.LengthMetric;
import org.biojava.bio.structure.align.benchmark.metrics.Metric;
import org.biojava.bio.structure.align.benchmark.metrics.PercentCorrectMetric;
import org.biojava.bio.structure.align.benchmark.metrics.RMSDMetric;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * An AlignmentStats object holds the result of calling the 
 * {@link Metric#calculate(MultipleAlignment, AFPChain, Atom[], Atom[]) calculate}
 * method for a series of {@link Metric}s on a particular pair of structures.
 * The result is cached and methods for outputting the stats are provided.
 * @author Spencer Bliven
 *
 */
public class AlignmentStats  {

	private List<Metric> metrics;
	private List<Double> results;
	private String name1, name2; //PDB identifiers
	
	public AlignmentStats(List<Metric> metrics, MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
		this.metrics = metrics;
		this.results = null;
		calculateStats(reference,align, ca1, ca2);
		
		String[] pdbIDs = reference.getNames();
		this.name1=pdbIDs[0];
		this.name2=pdbIDs[1];
	}
	
	/**
	 * @return the metrics
	 */
	public List<Metric> getMetrics() {
		return metrics;
	}

	public void writeHeader(Writer w) throws IOException {
		AlignmentStats.writeHeader(w, this.metrics);
	}
	
	public static void writeHeader(Writer w, List<Metric> metrics) throws IOException {
		w.write("PDB1\tPDB2");

		if(metrics.size()<1) {
			return;
		}
		
		//Print header
		for(int i=0;i<metrics.size();i++) {
			w.write('\t');
			w.write(metrics.get(i).getName());
		}
		w.write('\n');	
	}
	
	public void writeStats(Writer w) throws IOException {
		if(results == null) {
			throw new IllegalStateException("calculateStats() has not been run.");
		}
		
		assert( metrics.size() == results.size() );
		
		//Print identifiers
		w.write(name1);
		w.write('\t');
		w.write(name2);
		
		//Print each metric
		Iterator<Double> objIt = results.iterator();
		Iterator<Metric> metricIt = metrics.iterator();
		
		while(metricIt.hasNext()) {
			w.write('\t');
			w.write(metricIt.next().format(objIt.next()));
		}
		w.write('\n');
		assert(!objIt.hasNext());
	}
	
	protected void calculateStats(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2) {
		results = new ArrayList<Double>(metrics.size());
		
		for(Metric m : metrics) {
			results.add(m.calculate(reference, align, ca1, ca2));
		}
	}
	
	/**
	 * A fairly comprehensive list of metrics to compute
	 * @return
	 */
	public static List<Metric> defaultMetrics() {
		ArrayList<Metric> metrics = new ArrayList<Metric>();
		metrics.add(new LengthMetric.Reference());
		metrics.add(new LengthMetric.Alignment());
		metrics.add(new RMSDMetric.Reference());
		metrics.add(new RMSDMetric.Alignment());
		metrics.add(new PercentCorrectMetric());
		return metrics;
	}
}
