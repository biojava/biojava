package org.biojava3.genome.parsers.gff;

import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * A Feature corresponds to a single row in a GFF file.
 *
 * @author Hanno Hinsch
 */
public class Feature implements FeatureI {

    private Location mLocation;
    private String mSeqname;
    private String mSource;
    private String mType;
    private double mScore;			//or . if none
    private int mFrame;				//0,1,2
    private String mAttributes;			//any trailing stuff
    private HashMap<String, String> mUserMap;

    /**
     * Get the sequence name. (GFF field 1). Note that feature objects have
     * no link or reference to the actual sequence object to which
     * they refer; they are completely uncoupled.
     *
     * @return Sequence name.
     */
    public String seqname() {
        return mSeqname;
    }

    ;

    /**
     * Get source (aka method). (GFF field 2). This is often the name of
     * the program or procedure that created the features.
     *
     * @return Source field.
     */
    public String source() {
        return mSource;
    }

    ;

    /**
     * Get feature type, such as "exon" or "CDS". (GFF field 3).
     *
     * @return Feature type.
     */
    public String type() {
        return mType;
    }

    ;

    /**
     * Get location of feature. Note that feature objects have
     * no link or reference to the actual sequence object to which
     * they refer; they are completely uncoupled.
     *
     * @return Location of feature.
     */
    public Location location() {
        return mLocation;
    }

    /**
     * Get score. (GFF field 7). The meaning of the score varies from file to file.
     *
     * @return Score value. 
     */
    public double score() {
        return mScore;
    }

    ;

    /**
     * Get frame (aka phase). (GFF field 8). Specifies the offset of the
     * first nucleotide of the first in-frame codon, assuming this feature 
     * is a dna/rna sequence that codes
     * for a protein. If you
     * intend to use this field, you probably want to look it up on the web first.
     *
     * @return The frame (0, 1, 2).
     */
    public int frame() {
        return mFrame;
    }

    ;

    /**
     * Get the string of key/value attributes. (GFF field 9). The format and
     * meaning of this field varies from flavor to flavor of GFF/GTF. This method
     * simply returns the whole string. Other methods in this class make assumptions
     * about its format and provide additional utility.
     *
     * @return The attribute string.
     */
    public String attributes() {
        return mAttributes;
    }

    ;

    private Feature() {
    }

    ;        //unavailable

    /**
     * Make a copy of the specified feature. The mappings in the userMap() HashMap
     * are copied, so each feature has independent user data. Note, however, that the
     * actual objects in the HashMap are shared (not copied), so a change to such an object may
     * affect multiple features.
     *
     * @param feature Feature to clone.
     */
    public Feature(Feature feature) {

        mSeqname = feature.mSeqname;
        mSource = feature.mSource;
        mType = feature.mType;
        mLocation = feature.mLocation;
        mScore = feature.mScore;
        mFrame = feature.mFrame;
        mAttributes = feature.mAttributes;
        initAttributeHashMap();
        mUserMap = new HashMap<String, String>(feature.mUserMap);
    }

    /**
     * Construct a new Feature from raw data (usually a GFF row).
     *
     * @param seqname The sequence name field (field 1).
     * @param source The source or method field (field 2).
     * @param type The type of feature field (field 3).
     * @param location The location of the feature. (calculated from GFF start, end and strand fields).
     * @param score The score field (field 7).
     * @param frame The frame or phase field (field 8).
     * @param attributes A string of key/value pairs separated by semicolons (field 9).
     */
    public Feature(String seqname, String source, String type, Location location, Double score, int frame, String attributes) {

        mSeqname = seqname;
        mSource = source;
        mType = type;
        mLocation = location;
        mScore = score;
        mFrame = frame;
        mAttributes = attributes;
        initAttributeHashMap();
        mUserMap = new HashMap<String, String>();

    }

    /**
     * Get HashMap of user data. Each Feature object has a Java HashMap object
     * which can be used to annotate the Feature. JavaGene does not use or interpret
     * the keys or values. The values can be any subtype of the Java Object class.
     *<br><br>
     * If a Feature is constructed from data fields, the initial HashMap has no mappings (is empty).
     * If a Feature is constructed from another Feature, a copy of the mappings is made.
     * Note that the Objects in the copied mapping are shared, even though the mapping itself
     * is copied (not shared). Thus removing or adding a mapping to one Feature will not affect the
     * other, but changing an Object which is part of an established mapping may affect both Features.
     *
     * @return The user HashMap.
     */
    public HashMap<String, String> userData() {
        return mUserMap;
    }

     HashMap<String,String> attributeHashMap = new HashMap<String,String>();
    
    private void initAttributeHashMap(){
       String[] values = mAttributes.split(";");
       for(String attribute : values){
           attribute = attribute.trim();
           int equalindex = attribute.indexOf("=");
           String splitData = "=";
           if(equalindex == -1) //gtf uses space and gff3 uses =
               splitData = " ";
           String[] data = attribute.split(splitData);
           String value = "";
           if(data.length >= 2 && data[1].indexOf('"') != -1){ // an attibute field could be empty
               value = data[1].replaceAll('"' + "","").trim();
           }else if(data.length >= 2){
               value = data[1].trim();
           }
           attributeHashMap.put(data[0].trim(), value);
       }
    }
    
    /**
     * Get value of specified attribute key. Returns null if the attribute key has no value (does not exist).
     * Keys are case-sensitive. Assumes attributes are correctly formatted in GFF style.
     * Known bug: a semicolon within a quoted value will cause parse failure.
     *
     * @param key The key.
     * @return The corresponding value. Null if the key has no value defined.
     */
    public String getAttribute(String key) {
        
        return attributeHashMap.get(key);
    }

    public String getAttributeOld(String key) {
        int start = 0;

        int end = mAttributes.indexOf(';');
        while (0 < end) {
            //find the first word (up to space) in chunk,
            // see if it is this key
            int i = mAttributes.indexOf(' ', start);
            if (0 < i && i < end) {
                if (mAttributes.substring(start, i).equals(key)) {
                    //remove quotes, if needed
                    if (mAttributes.charAt(i + 1) == '\"' && mAttributes.charAt(end - 1) == '\"') {
                        return mAttributes.substring(i + 2, end - 1);//return attribute
                    } else {
                        return mAttributes.substring(i + 1, end);	//return attribute
                    }
                }
            }
            start = end + 2;	//skip required semicolon and single space
            end = mAttributes.indexOf(';', start);
        }

        return null;
    }

    public boolean hasAttribute(String key) {
        return attributeHashMap.containsKey(key);
    }

    public boolean hasAttribute(String key, String value) {
        String data = getAttribute(key);
        if(data == null)
            return false;
        if(data.equals(value))
            return true;
        else
            return false;
    }

    /**
     * Get the first item (everything before first semicolon, if it has one)
     * in the attribute field, which is assumed to
     * be a group identifer. This is appropriate for GFF1 files and variants. It is not
     * appropriate for GTF and GFF2 files, although they may use a named attribute key,
     * such as "gene_id" or "transcript_id", for grouping.
     *
     * @return The group id. Everything before the first semicolon in the attributes string (minus trailing whitespace).
     */
    public String group() {
        int i = mAttributes.indexOf(';');
        return (i < 0) ? mAttributes.trim() : mAttributes.substring(0, i).trim();
    }

    /**
     *
     */
    public String toString() {
        String s = mSeqname + '\t';
        s += mSource + '\t';
        s += mType + '\t';
        s += mLocation.start() + "\t";
        s += mLocation.end() + "\t";
        s += Double.toString(mScore) + "\t";

        if (mFrame == -1) {
            s += ".\t";
        } else {
            s += mFrame + "\t";
        }

        s += mAttributes;

        return s;
    }

    /**
     * @deprecated
     */
    public static void main(String args[])
            throws Exception {
        //Feature f= new Feature();
        //intentionally perverse
        //f.group= "gene_id transcript; transcript \"gene_id fantom2\"; ";
        //	f.addAttribute( "author", "julian" );
        //	f.addAttribute( "curator", "nick" );
        //	f.addAttribute( "author", "hanno" );
        //Log.log( f.group );
        //f.addAttribute( "perverse", "foo;goo" );
        //assert f.getAttribute( "perverse").equals( "foo;goo" );
        //	assert f.getAttribute( "gene_id" ).equals( "transcript" );
        //	assert f.getAttribute( "author" ).equals( "julian hanno" );
        //	assert f.getAttribute( "curator" ).equals( "nick" );
        //	assert f.getAttribute( "transcript").equals( "gene_id fantom2" );
        //Log.log( "passed test." );
    }

	@Override
	public HashMap<String, String> getAttributes() {
		
		return attributeHashMap;
	}
}
