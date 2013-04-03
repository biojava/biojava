/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.benchmark.metrics.ConsistencyMetric;
import org.biojava.bio.structure.align.benchmark.metrics.AlignmentLengthMetric;
import org.biojava.bio.structure.align.benchmark.metrics.MetaDataMetric;
import org.biojava.bio.structure.align.benchmark.metrics.Metric;
import org.biojava.bio.structure.align.benchmark.metrics.PercentCorrectMetric;
import org.biojava.bio.structure.align.benchmark.metrics.ProteinLengthMetric;
import org.biojava.bio.structure.align.benchmark.metrics.RMSDMetric;
import org.biojava.bio.structure.align.benchmark.metrics.TMScoreMetric;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * An AlignmentStats object holds the result of calling the 
 * {@link Metric#calculate(MultipleAlignment, AFPChain, Atom[], Atom[], Map) calculate}
 * method for a series of {@link Metric}s on a particular pair of structures.
 * The result is cached and methods for outputting the stats are provided.
 * @author Spencer Bliven
 *
 */
public class AlignmentStats  {

	private List<Metric> metrics;
	private List<Double> results;
	private String name1, name2; //PDB identifiers
	
	/**
	 * Create a new AlignmentStats
	 * @param metrics The set of metrics to evaluate for the aligment pair
	 * @param reference The "reference" alignment (must be pairwise)
	 * @param align The other alignment, which will be compared to the reference (must be pairwise)
	 * @param ca1 All CA atoms of protein 1
	 * @param ca2 All CA atoms of protein 2
	 * @param metaData A map storing arbitrary objects relating to the alignment process. May be null.
	 * The {@link MetaDataMetric} displays values from here. Currently used to store "alignmentTime".
	 */
	public AlignmentStats(List<Metric> metrics, MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		this.metrics = metrics;
		this.results = null;
		calculateStats(reference,align, ca1, ca2, metaData);
		
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

	public double getResult(int i) {
		if(results == null) {
			throw new IllegalStateException("calculateStats() has not been run.");
		}
		return results.get(i);
	}
	
	/**
	 * @return the name1
	 */
	public String getName1() {
		return name1;
	}

	/**
	 * @return the name2
	 */
	public String getName2() {
		return name2;
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
	
	protected void calculateStats(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String,Object> metaData) {
		results = new ArrayList<Double>(metrics.size());
		
		for(Metric m : metrics) {
			results.add(m.calculate(reference, align, ca1, ca2, metaData));
		}
	}
	
	/**
	 * A fairly comprehensive list of metrics to compute
	 * @return
	 */
	public static List<Metric> defaultMetrics() {
		ArrayList<Metric> metrics = new ArrayList<Metric>();
		metrics.add(new ProteinLengthMetric(0));
		metrics.add(new ProteinLengthMetric(1));
		metrics.add(new MetaDataMetric("alignmentTime","time"));
		metrics.add(new AlignmentLengthMetric.Reference());
		metrics.add(new AlignmentLengthMetric.Alignment());
		metrics.add(new RMSDMetric.Reference());
		metrics.add(new RMSDMetric.Alignment());
		metrics.add(new TMScoreMetric.Reference());
		metrics.add(new TMScoreMetric.Alignment());
		//metrics.add(new PercentCorrectMetric());
		metrics.add(new ConsistencyMetric() );
		metrics.add(new ConsistencyMetric(4) );
		return metrics;
	}
}
