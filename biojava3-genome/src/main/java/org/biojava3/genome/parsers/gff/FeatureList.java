package org.biojava3.genome.parsers.gff;

import java.util.*;
import java.util.Map.Entry;

import org.biojava3.core.sequence.DNASequence;



/**
 * A list of FeatureI objects implemented using a Java ArrayList; corresponds to a GFF file.
 * This class is implemented entirely using FeatureI objects, so everything here will work
 * correctly if you choose to implement your own feature class -- there are no dependencies
 * on JavaGene's native Feature class.
 *
 *
 * @author Hanno Hinsch, Carmelo Foti
 */
public class FeatureList extends ArrayList<FeatureI> {

	 Map<String, Map<String,List<FeatureI>>> featindex = new HashMap<String,Map<String,List<FeatureI>>>();
	Location mLocation;			//genomic location (union of feature locations)

	/**
	 * Construct an empty list.
	 */
	public FeatureList() {
		mLocation = null;
	}

	/**
	 * Construct a new list containing the same features
	 * as the specified list.
	 *
	 * @param features An existing list or collection of FeatureI objects.
	 */
	public FeatureList(Collection<FeatureI> features) {
		this();
		this.add(features);

		mLocation = null;
	}

	/**
	 * Add specified feature to the end of the list. Updates the bounding location of the
	 * feature list, if needed.
	 *
	 * @param feature The FeatureI object to add.
	 * @return True if the feature was added.
	 */
	public boolean add(FeatureI feature) {
		if (mLocation == null) {
			mLocation = feature.location().plus();
		} else if (null != feature.location()) {
			mLocation = mLocation.union(feature.location().plus());
		}
		for (Entry<String, String> entry : feature.getAttributes().entrySet()){
			if (featindex.containsKey(entry.getKey())){
				Map<String,List<FeatureI>> feat = featindex.get(entry.getKey());
				if (feat==null){
					feat= new HashMap<String,List<FeatureI>>();
				}
				List<FeatureI> features = feat.get(entry.getValue());
				if (features==null){
					features = new ArrayList<FeatureI>();
				}
				features.add(feature);
				feat.put(entry.getValue(), features);
				featindex.put(entry.getKey(), feat);
				//featindex.put(key, value)
			}
		}

		return super.add(feature);
	}

	/**
	 * Add all features in the specified list or collection to this list.
	 *
	 * @param list The collection of FeatureI objects.
	 */
	public void add(Collection<FeatureI> list) {
		for (FeatureI f : list) {
			add(f);
		}
	}

	/**
	 * The union of all locations of all features in this list, mapped to the positive strand.
	 * If an added feature is on the negative strand, its positive strand image is added
	 * to the union.
	 * The bounding location is not updated when a feature is removed from the list, so
	 * it is not guaranteed to be the minimal bounding location.
	 *
	 * @return A location that is the union of all feature locations in the list.
	 */
	public Location bounds() {
		return mLocation;
	}

	/**
	 * Check size of gaps between successive features in list. The features in
	 * the list are assumed to be appropriately ordered.
	 *
	 * @param gapLength The minimum gap length to consider. Use a gapLength
	 * of 0 to check if features are contiguous.
	 * @return True if list has any gaps equal to or greater than gapLength.
	 */
	public boolean hasGaps(int gapLength) {
		Location last = null;
		for (FeatureI f : this) {
			if (last != null && gapLength <= f.location().distance(last)) {
				return true;
			} else {
				last = f.location();
			}
		}

		return false;
	}




	/**
	 * Concatenate successive portions of the specified sequence
	 * using the feature locations in the list. The list is assumed to be appropriately
	 * ordered.
	 *
	 * @param sequence The source sequence from which portions should be selected.
	 * @return The spliced data.
	 * @throws IllegalStateException Out of order or overlapping FeatureI locations detected.
	 *
	 */
	public String splice(DNASequence sequence) {
		StringBuilder subData = new StringBuilder();		
		Location last = null;

		for (FeatureI f : this) {
			Location loc = f.location();

			if (last == null || loc.startsAfter(last)) {
				subData.append(sequence.getSubSequence(loc.start(), loc.end()).toString());
				last = loc;
			} else {
				throw new IllegalStateException("Splice: Feature locations should not overlap.");
			}

		}

		return subData.toString();
	}

	/**
	 * Create a collection of all unique group ids in the list, as defined
	 * by the group() method of the features. For example, if the
	 * features are from a GFF1 file, then each group id identifies a particular gene,
	 * and this method returns a collection of all gene ids.
	 *
	 * @return A collection (suitable for iteration using Java's "for" loop) of all the
	 * group ids found in this list. The order of the values is undefined; it will not match
	 * the order of features in the list.
	 */
	public Collection<String> groupValues() {
		Set<String> set = new HashSet<String>();
		for (FeatureI f : this) {
			//enter in a set -- removes duplicates
			set.add(f.group());
		}

		return set;
	}

	/**
	 * Create a collection of the unique values for the specified key.
	 * Example: For GTF files, using the "gene_id" key will give the names of all
	 * the genes in this list.
	 *
	 * @return A collection (suitable for iteration using Java's "for" loop) of all the
	 * values found for this key. The order of the values is undefined; it will not match
	 * the order of features in the list.
	 */
	public Collection<String> attributeValues(String key) {
		if (featindex.containsKey(key)){
			Map<String, List<FeatureI>> map = featindex.get(key);
			Collection<String> result = map.keySet();
			if (result == null) result = new HashSet<String>();
			return Collections.unmodifiableCollection(result);
		}
		LinkedHashMap<String, String> hash = new LinkedHashMap<String, String>();
		for (FeatureI f : this) {
			//enter as a key -- removes duplicates
			hash.put(f.getAttribute(key), null);
		}

		return Collections.unmodifiableCollection(hash.keySet());
	}

	/**
	 * Create a list of all features that have the specified group id, as defined by
	 * the group() method of the features. 
	 *
	 * @param groupid The group to match.
	 * @return A list of features having the specified group id.
	 */
	public FeatureList selectByGroup(String groupid) {
		FeatureList list = new FeatureList();
		for (FeatureI f : this) {
			if (f.group().equals(groupid)) {
				list.add(f);
			}
		}

		return list;
	}

	/**
	 * Create a list of all features that are of the specified type, as defined by
	 * the type() method of the features. 
	 * This might be, for example, "exon" or "CDS".
	 *
	 * @param type The type to match.
	 * @return A list of features of the specified type.
	 */
	public FeatureList selectByType(String type) {
		FeatureList list = new FeatureList();
		for (FeatureI f : this) {
			if (f.type().equals(type)) {
				list.add(f);
			}
		}

		return list;
	}

	/**
	 * Create a list of all features that include the specified attribute key/value pair.
	 *
	 * @param key The key to consider.
	 * @param value The value to consider.
	 * @return A list of features that include the key/value pair. 
	 */
	public FeatureList selectByAttribute(String key, String value) {
		if (featindex.containsKey(key)){
			Map<String,List<FeatureI>> featuresmap = featindex.get(key);
			if (featuresmap==null) return new FeatureList();
			List<FeatureI> list = featuresmap.get(value);
			if (list == null){
				return new FeatureList();
			}
			return  new FeatureList(list);
		}
		FeatureList list = new FeatureList();
		for (FeatureI f : this) {
			if (f.hasAttribute(key, value)) {
				list.add(f);
			}
		}
		return list;
	}
	/**
	 * Create a list of all features that include the specified attribute key.
	 *
	 * @param key The key to consider.
	 * @return A list of features that include the key. 
	 */
	public FeatureList selectByAttribute(String key) {
		FeatureList list = new FeatureList();
		if (featindex.containsKey(key)){
			Map<String, List<FeatureI>> featsmap =featindex.get(key);
			for (List<FeatureI> feats: featsmap.values()){
				list.addAll(Collections.unmodifiableCollection(feats));
			}
			return list;
		}
		
		for (FeatureI f : this) {
			if (f.hasAttribute(key)) {
				list.add(f);
			}
		}
		return list;
	}

	/**
	 * Create a list of all features that include the specified key/value pair in their userMap().
	 *
	 * @param key The key to consider.
	 * @param value The value to consider.
	 * @return A list of features that include the key/value pair. 
	 */
	public FeatureList selectByUserData(String key, Object value) {
		FeatureList list = new FeatureList();
		for (FeatureI f : this) {
			Object o = f.userData().get(key);
			if (o != null && o.equals(value)) {
				list.add(f);
			}
		}
		return list;
	}
	/**
	 * Create a list of all features that include the specified key in their userMap().
	 *
	 * @param key The key to consider.
	 * @return A list of features that include the key. 
	 */
	public FeatureList selectByUserData(String key) {
		FeatureList list = new FeatureList();
		for (FeatureI f : this) {
			if (f.userData().containsKey(key)) {
				list.add(f);
			}
		}
		return list;
	}

	/**
	 * Create a list of all features that overlap the specified location on the specified
	 * sequence.
	 *
	 * @param seqname The sequence name. Only features with this sequence name will be checked for overlap.
	 * @param location The location to check.
	 * @param useBothStrands If true, locations are mapped to their positive strand image
	 * before being checked for overlap. If false, only features whose locations are 
	 * on the same strand as the specified location will be considered for inclusion.
	 * @return The new list of features that overlap the location.
	 */
	public FeatureList selectOverlapping(String seqname, Location location, boolean useBothStrands)
			throws Exception {
		FeatureList list = new FeatureList();

		for (FeatureI feature : this) {
			boolean overlaps = false;
			if (feature.seqname().equals(seqname)) {
				if (location.isSameStrand(feature.location())) {
					overlaps = feature.location().overlaps(location);
				} else if (useBothStrands) {
					overlaps = feature.location().overlaps(location.opposite());
				}
			}
			if (overlaps) {
				list.add(feature);
			}
		}
		return list;
	}

	/**
	 * Create a list of all features that do not overlap the specified location on the specified sequence.
	 *
	 * @param seqname The sequence name. Only features with this sequence name will be checked for overlap.
	 * @param location The location to check.
	 * @param useBothStrands If true, locations are mapped to their positive strand image
	 * before being checked for overlap. If false, all features whose locations are
	 * on the opposite strand from the specified location will be considered non-overlapping.
	 * @return The new list of features that do not overlap the location.
	 */
	public FeatureList omitOverlapping(String seqname, Location location, boolean useBothStrands) {
		FeatureList list = new FeatureList();

		for (FeatureI feature : this) {
			boolean overlaps = false;
			if (feature.seqname().equals(seqname)) {
				if (location.isSameStrand(feature.location())) {
					overlaps = feature.location().overlaps(location);
				} else if (useBothStrands) {
					overlaps = feature.location().overlaps(location.opposite());
				}
			}

			if (!overlaps) {
				list.add(feature);
			}
		}

		return list;
	}

	/**
	 * Check if any feature in list has the specified attribute key.
	 *
	 * @param key The attribute key to consider.
	 * @return True if at least one feature has the attribute key.
	 */
	public boolean hasAttribute(String key) {
		if (featindex.containsKey(key)){
			Map<String, List<FeatureI>> mappa = featindex.get(key);
			if (mappa!= null && mappa.size()>0)return true;
			return false;
		}
		for (FeatureI f : this) {
			if (f.hasAttribute(key)) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Check if any feature in list has the specified attribute key/value pair.
	 *
	 * @param key The attribute key to consider.
	 * @param value The attribute value to consider.
	 * @return True if at least one feature has the key/value pair.
	 */
	public boolean hasAttribute(String key, String value) {
		if (featindex.containsKey(key)){
			Map<String, List<FeatureI>> mappa = featindex.get(key);
			if (mappa == null) return false;
			if (mappa.containsKey(value)) return true;
			return false;
		}

		for (FeatureI f : this) {
			if (f.hasAttribute(key, value)) {
				return true;
			}
		}

		return false;
	}

	/**
	 * Return a string representation of all features in this list.
	 *
	 * @return A string.
	 */
	public String toString() {
		StringBuilder  s = new StringBuilder("FeatureList: >>\n");
		for (FeatureI f : this) {
			s.append( f.seqname() + ":" + f.toString() + "\n");
		}

		s.append("\n<<\n");
		return s.toString();
	}

	/**
	 * used by sort routine
	 */
	private class FeatureComparator implements Comparator<FeatureI> {

		public int compare(FeatureI a, FeatureI b) {
			if (a.seqname().equals(b.seqname()) && a.location().isSameStrand(b.location())) {
				return a.location().start() - b.location().start();		//sort on start
			} else {
				throw new IndexOutOfBoundsException("Cannot compare/sort features whose locations are on opposite strands or with different seqname().\r\n" + a.toString() + "\r\n" + b.toString() );
			}
		}
	}

	/**
	 * Create a new list that is ordered by the starting index of the features' locations. All locations
	 * must be on the same strand of the same sequence.
	 *
	 * @return An ordered list.
	 * @throws IndexOutOfBoundsException Cannot compare/sort features whose locations are on opposite strands, or
	 * whose seqnames differ.
	 */
	public FeatureList sortByStart() {
		FeatureI array[] = (FeatureI[]) toArray(new FeatureI[1]);

		Arrays.sort(array, new FeatureComparator());

		return new FeatureList(Arrays.asList(array));
	}

	/**
	 * @deprecated
	 *
	 */
	// FIXME features may have a null location() !!
	static public void main(String args[]) {
	}

	
	/**
	 * Add a list of attributes that will be used as indexes for queries 
	 * @param indexes  the List containing the attribute_id
	 */
	public void addIndexes(List<String> indexes) {
		for (String index : indexes){
			addIndex(index);
		}

	}
	/**
	 * Add an attribute that will be used as index for queries 
	 * @param index an attribute_id
	 */
	public void addIndex(String index) {
		featindex.put(index, null);
	}
}
