/*
 *                  BioJava development code
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
 * Created on Dec 5, 2007
 * 
 */

package org.biojava.dasobert.feature;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;



/** converts the features from their "raw" representation as a Map
 * into a Feature class
 * @author Andreas Prlic
 *
 */
public class FeatureTrackConverter {

	public static final Color HELIX_COLOR  = new Color(255,51,51);
	public static final Color STRAND_COLOR = new Color(255,204,51);
	public static final Color TURN_COLOR   = new Color(204,204,204); 


	// some  annotation types, for which there is a special treatment
	public static final String DISULFID_TYPE = "DISULFID";
	public static final String SECSTRUC_TYPE = "SECSTRUC";
	public static final String METAL_TYPE    = "METAL";
	public static final String MSD_SITE_TYPE = "MSD_SITE";

	String type;
	public static final String TYPE_HISTOGRAM = "histogram";
	public static final String TYPE_DEFAULT   = "default";


	// rotate between these colors
	public static final  Color[] entColors = new Color []{
		new Color(51,51,255), // blue
		new Color(102,255,255),    // cyan
		new Color(153,255,153), // green
		new Color(153,255,153), // green
		new Color(255,153,153), // pink
		new Color(255,51,51),   // red
		new Color(255,51,255)    // pink 
	};

	public static final String[] txtColors = new String[] { 
		"blue",
		"pink",
		"green",
		"yellow",
		"red",
		"cyan",
	"pink"};

	Map[] stylesheet;

	boolean isHistogram = false;


	public FeatureTrackConverter(){
		type = TYPE_DEFAULT;	
		stylesheet = new Map[0];
	}

	public FeatureTrackConverter(Map[] stylesheet){

		if (stylesheet == null)
			stylesheet = new Map[0];

		this.stylesheet = stylesheet;

	}


	public FeatureTrackConverter(Map[] stylesheet, boolean isHistogram){
		this(stylesheet);
		this.isHistogram = isHistogram;
		if ( isHistogram)
			type = TYPE_HISTOGRAM;
	}
	public String getType() {
		return type;
	}

	public void setType(String type) {
		if ( type.equals(TYPE_HISTOGRAM))
			isHistogram = true;
		this.type = type;
	}

	public boolean isHistogram() {
		return isHistogram;
	}


	public void setHistogram(boolean isHistogram) {
		this.isHistogram = isHistogram;
	}

	

	

	public  FeatureTrack[] convertMap2Features(Map<String,String>[] mapfeatures){
		List<FeatureTrack> features = new ArrayList<FeatureTrack>();

		boolean first = true ;
		boolean secstruc = false ;
		boolean isGroup  = false;

		FeatureTrack feat = null ;
		Segment segment   = null ;

		int featuresCounter = 0;
		String prevGroup = null;

		for (int i = 0 ; i< mapfeatures.length;i++) {

			Map<String,String> currentFeatureMap = mapfeatures[i];
			String type = (String) currentFeatureMap.get("TYPE") ;

			
			
			String group = currentFeatureMap.get("GROUP");
			if ( group != null){
				if ( prevGroup != null) {
					if ( group.equals(prevGroup)){
						feat.setName(group);
						isGroup = true;
					} else {
						isGroup = false;
					}
				} else {
					isGroup = false;
				}
			} else {
				isGroup = false;
			}			

			// we are skipping literature references for the moment 
			// TODO: add a display to spice for non-positional features
			//
			if ( type.equals("reference") || type.equals("GOA")){
				continue ;
			}

			if (! first) 
			{
				// if not first feature

				if ( (! secstruc ) && (! isGroup) )  {             

					// if not secondary structure and not in a group ...
					features = testAddFeatures(features,feat);


				} else if ( ! 
						(
								type.equals("HELIX")  || 
								type.equals("STRAND") || 
								type.equals("TURN")  
						) 
				)
				{
					// end of secondary structure
					secstruc = false ;
					if ( feat != null && (! isGroup)) {						
						features = testAddFeatures(features,feat);                        
					}

				}
			} // end of not first

			first = false ;             
			if ( (! secstruc) && (! isGroup)) {				
				featuresCounter +=1;
				feat = getNewFeat(currentFeatureMap);       
			}




			if (type.equals("STRAND")){
				secstruc = true ;				
				currentFeatureMap.put("colorTxt","yellow");
				feat.setName("SECSTRUC");       
				feat.setType("SECSTRUC");
			}

			else if (type.equals("HELIX")) {
				secstruc = true ;				
				currentFeatureMap.put("colorTxt","red");
				feat.setName("SECSTRUC");
				feat.setType("SECSTRUC");
			}   

			else if (type.equals("TURN")) {
				secstruc = true ;				
				currentFeatureMap.put("colorTxt","white");

				feat.setName("SECSTRUC");
				feat.setType("SECSTRUC");
			}     
			else {
				secstruc = false ;				
				currentFeatureMap.put("colorTxt",txtColors[featuresCounter%txtColors.length]);
				if ( ! isGroup) {
					try {
						feat.setName(type);
					
					} catch ( NullPointerException e) {
						//e.printStackTrace();
						feat.setName("null");
					}
				}
			}

			segment = getNewSegment(currentFeatureMap);

			feat.addSegment(segment);       
			prevGroup = group;
		}   

		if ( feat != null)  
			features =testAddFeatures(features,feat);


		return (FeatureTrack[]) features.toArray(new FeatureTrack[features.size()]);
	}

	/**	 test if this features is added as a new feature to the features list, 
	 * or if it is joint with an already existing one...
	 * 
	 * @param features
	 * @param newFeature
	 * @return a List of FeatureTrack objects
	 */
	protected  List<FeatureTrack> testAddFeatures(List<FeatureTrack> features,FeatureTrack newFeature){
		
		//System.out.println("testing " + newFeature + " " + newFeature.getScore());   
		Iterator iter = features.iterator();


		if ( isHistogramFeatureType(newFeature)) {  
			
			// return histogram type features
			type = TYPE_HISTOGRAM;

			Segment seg = getHistogramSegmentFromFeature(newFeature);

			while (iter.hasNext()){
				FeatureTrack knownFeature = (FeatureTrack) iter.next() ;
				String knownType = knownFeature.getType();

				//System.out.println("found histogram style " + feat);
				// set type of this DAS source to being HISTOGRAM style


				if ( knownType.equals(newFeature.getType())){
					// convert the feature into a HistogramSegment and add to the already known feature

					knownFeature.addSegment(seg);
					// we can return now
					return features;
				}


			}
			// we could not link this to any existing feature
			// convert it to a new HistogramFeature
			HistogramFeature hfeat = new HistogramFeature();

			hfeat.setLink(newFeature.getLink());
			hfeat.setMethod(newFeature.getMethod());
			hfeat.setName(newFeature.getName());
			hfeat.setNote(newFeature.getNote());
			hfeat.setScore("0");
			hfeat.setSource(newFeature.getSource());
			hfeat.addSegment(seg);
			hfeat.setType(newFeature.getType());

			newFeature = hfeat;
			features.add(newFeature);
			return features;
		} 



		while (iter.hasNext()){
			FeatureTrack knownFeature = (FeatureTrack) iter.next() ;
			// this only compares method source and type ...
			boolean sameFeat = false;
			if ( knownFeature.equals(newFeature))
				sameFeat = true;

			if ( ( knownFeature.getSource().equals(newFeature.getSource() )) &&
					( knownFeature.getMethod().equals(newFeature.getMethod())) &&
					( knownFeature.getNote().equals(newFeature.getNote())) &&
					isSecondaryStructureFeat(knownFeature) && 
					isSecondaryStructureFeat(newFeature))
				sameFeat =true;

			if ( sameFeat) {

				// seems to be of same type, method and source, so check if the segments can be joined

				List<Segment> tmpsegs = knownFeature.getSegments();
				Iterator segiter = tmpsegs.iterator();
				List<Segment> newsegs = newFeature.getSegments();
				Iterator newsegsiter = newsegs.iterator();
				boolean overlap = false;
				while (newsegsiter.hasNext()){
					Segment newseg = (Segment)newsegsiter.next();


					while (segiter.hasNext()){
						Segment tmpseg = (Segment) segiter.next();

						if (  tmpseg.overlaps(newseg))
							overlap = true;
					}
				}

				if ( ! overlap){
					// add all new segments to old features...
					newsegsiter = newsegs.iterator();
					while (newsegsiter.hasNext()){
						Segment newseg = (Segment)newsegsiter.next();
						knownFeature.addSegment(newseg);
					}

					return features;
				} 
			}

		}

		//      if we get here, the  features could not be joint with any other one, so there is always some overlap
		// add to the list of known features
		features.add(newFeature);
		return features;
	}

	private FeatureTrack getNewFeat(Map currentFeatureMap) {
		FeatureTrack feat = new FeatureTrackImpl();
		//logger.finest(currentFeatureMap);
		//System.out.println("DrawableDasSource " + currentFeatureMap);
		feat.setSource((String)currentFeatureMap.get("dassource"));
		feat.setName(  (String)currentFeatureMap.get("NAME"));
		feat.setType(  (String)currentFeatureMap.get("TYPE"));
		feat.setLink(  (String)currentFeatureMap.get("LINK"));
		feat.setNote(  (String)currentFeatureMap.get("NOTE"));
		
		String typeID       = (String) currentFeatureMap.get("TYPE_ID");
		String typeCategory = (String) currentFeatureMap.get("TYPE_CATEGORY");
		feat.setTypeID(typeID);
		feat.setTypeCategory(typeCategory);
		
		String method = (String)currentFeatureMap.get("METHOD");
		if ( method == null) { method = "";}
		feat.setMethod(method);
		feat.setScore( (String)currentFeatureMap.get("SCORE"));
		return feat ;
	}

	private Segment getNewSegment(Map featureMap) {
		Segment s = new SegmentImpl();
		String sstart = (String)featureMap.get("START") ;
		String send   = (String)featureMap.get("END")   ;
		int start = Integer.parseInt(sstart) ;
		int end   = Integer.parseInt(send)   ;
		s.setStart(start);
		s.setEnd(end);
		s.setName((String)featureMap.get("TYPE"));
		s.setTxtColor((String)featureMap.get("colorTxt"));  
		s.setColor((Color)featureMap.get("color"));
		s.setNote((String) featureMap.get("NOTE"));
		return s ;

	}

	private boolean isSecondaryStructureFeat(FeatureTrack feat){
		String type = feat.getType();
		if (
				type.equals("HELIX")  || 
				type.equals("STRAND") || 
				type.equals("TURN")
		) return true;
		return false;
	}

	private boolean isHistogramFeatureType(FeatureTrack feat){
		String ftype = feat.getType();

		Map[] style = stylesheet;

		//System.out.println("is HistogramFeature type " + ftype + " " + style );


		// todo : move this info into a config file...

		if ( ftype.equals("hydrophobicity")){
			return true;
		}
		if ( getType().equals(TYPE_HISTOGRAM) )
			return true;



		if (style != null ) {

			for ( int i =0; i< style.length ; i++){
				Map m = style[i];

				// make sure the stylesheet is for this feature type
				String styleType = (String) m.get("type");
				if ( styleType != null) {
					if ( ! styleType.equals(ftype)){
						continue;
					}
				} else {
					continue;
				}

				String type = (String) m.get("style");
				if ( type != null) {
					//System.out.println("stylesheet type " + type);
					if ( type.equals("gradient") || ( type.equals("lineplot")) || ( type.equals("histogram"))){

						return true;
					}
				}
			}
		}

		return false;
	}


	private HistogramSegment getHistogramSegmentFromFeature(FeatureTrack feat){
		HistogramSegment s = new HistogramSegment();

		double score = 0.0;

		try {
			score = Double.parseDouble(feat.getScore());

		} catch (Exception e){
			//e.printStackTrace();
		}
		s.setScore(score);		
		List segments = feat.getSegments();
		if (segments.size() > 0){
			Segment seg = (Segment) segments.get(0);
			s.setName(seg.getName());
			s.setStart(seg.getStart());
			s.setEnd(seg.getEnd());
			s.setNote(seg.getNote());
			s.setColor(seg.getColor());
			s.setTxtColor(seg.getTxtColor());
		}


		return s;
	}


}
