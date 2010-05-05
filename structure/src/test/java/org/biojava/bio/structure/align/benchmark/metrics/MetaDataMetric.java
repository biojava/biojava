/**
 * 
 */
package org.biojava.bio.structure.align.benchmark.metrics;

import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.benchmark.MultipleAlignment;
import org.biojava.bio.structure.align.model.AFPChain;

/**
 * This Metric makes use of the metaData map passed to {@link #calculate()}
 * from the alignment step.
 * It allows the display of any key in the map.
 * <p>
 * Keys may be {@link StructureAlignment}-specific. Current general keys are:
 * <ul><li>alignmentTime</li>
 * </ul>
 * <p>
 * This class is a bit complex because it needs to be able to display arbitrary objects.
 * See {@link #setFormatter(MetaDataFormatter)} and {@link #setFormat(String)}.
 * 
 * @author Spencer Bliven
 *
 */
public class MetaDataMetric extends Metric {
	private String key;
	private String name;
	private String format;
	private MetaDataFormatter formatter;

	/**
	 * Create a MetaDataMetric with the specified key.
	 * Use the key as the display name.
	 */
	public MetaDataMetric(String key) {
		this(key,key);
	}

	/**
	 * Create a MetaDataMetric with the specified key and display name.
	 * Assumes that values for the key are Doubles.
	 * @param key
	 * @param name
	 */
	public MetaDataMetric(String key, String name) {
		this(key,name,new DoubleFormatter());
	}
	public MetaDataMetric(String key, String name, MetaDataFormatter formatter) {
		super();
		this.key = key;
		this.name = name;
		this.formatter = formatter;
		this.format = "%.04f";
	}



	/**
	 * @see org.biojava.bio.structure.align.benchmark.metrics.Metric#calculate(org.biojava.bio.structure.align.benchmark.MultipleAlignment, org.biojava.bio.structure.align.model.AFPChain, org.biojava.bio.structure.Atom[], org.biojava.bio.structure.Atom[], java.util.Map)
	 */
	@Override
	public double calculate(MultipleAlignment reference, AFPChain align,
			Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		if(metaData != null && metaData.containsKey(key)) {
			return formatter.format(metaData.get(key));
		} else {
			return 0.0;
		}
	}

	/**
	 * @see org.biojava.bio.structure.align.benchmark.metrics.Metric#getName()
	 */
	@Override
	public String getName() {
		return name;
	}
	
	@Override
	public String format(MultipleAlignment reference, AFPChain align, Atom[] ca1, Atom[] ca2, Map<String, Object> metaData) {
		if(metaData == null || ! metaData.containsKey(key)) {
			return "";
		} else {
			return this.format(calculate(reference, align, ca1, ca2, metaData));
		}
	}
	
	@Override
	public String format(double result) {
		return String.format(format,result);
	}

	/**
	 * @param name The description of this metric to use in table headers
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the key from metaData to display
	 */
	public String getKey() {
		return key;
	}

	/**
	 * @param key The key in the metaData to display
	 */
	public void setKey(String key) {
		this.key = key;
	}

	/**
	 * By default, MetaDataMetric just calls toString() on the object it gets
	 * from metaData. With this you can set a prettier format, such as "%02.4f" for doubles.
	 * 
	 * //@note This would be a lambda function in any other language in the world.
	 * @param format the format to set
	 * @see java.util.Formatter
	 */
	public void setFormat(String format) {
		this.format = format;
	}
	

	/**
	 * Set a custom formatter to handle a metaData value other than a Double.
	 * @param formatter the formatter to set
	 */
	public void setFormatter(MetaDataFormatter formatter) {
		this.formatter = formatter;
	}


	/**
	 * The metaData may contain objects of any type. These are converted into
	 * doubles by a MetaDataFormatter. Most of the time you should just store
	 * a double in the metaData map and use the default DoubleFormatter
	 * @author Spencer Bliven
	 *
	 */
	public static interface MetaDataFormatter {
		/**
		 * Formats the specified data
		 * @param data
		 * @throw ClassCastException if you use the wrong formatter for your metaData
		 * @return
		 */
		public double format(Object data);
	}
	/**
	 * A simple formatter that unwraps Doubles into doubles.
	 * @author Spencer Bliven
	 */
	public static class DoubleFormatter implements MetaDataFormatter {

		/**
		 * @see MetaDataFormatter#format(Object)
		 */
		//@Override
		public double format(Object data) {
			return ((Double) data).doubleValue();
		}
	}

	

}
