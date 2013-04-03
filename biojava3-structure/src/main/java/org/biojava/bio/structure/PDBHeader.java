package org.biojava.bio.structure;

import java.io.Serializable;
import java.lang.reflect.Method;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.biojava.bio.structure.quaternary.ModelTransformationMatrix;


/** A class that contains PDB Header information.
 *
 * @author Andreas Prlic
 * @since 1.6
 *
 */
public class PDBHeader implements PDBRecord, Serializable{
	/**
	 *
	 */
	private static final long serialVersionUID = -5834326174085429508L;
	String method;
	String title;
	String description;
	String idCode;
	String classification;
	Date depDate;
	Date modDate;
	String technique;
	float resolution;
	String authors;
	int nrBioAssemblies ;
	public static final float DEFAULT_RESOLUTION = 99;

	private Long id;
	public static final String newline = System.getProperty("line.separator");

	DateFormat dateFormat;

	Map<Integer,List<ModelTransformationMatrix>> tranformationMap ;

	public PDBHeader(){

		depDate = new Date(0);
		modDate = new Date(0);
		dateFormat = new SimpleDateFormat("dd-MMM-yy",Locale.US);
		resolution = DEFAULT_RESOLUTION;
		tranformationMap = new HashMap<Integer, List<ModelTransformationMatrix>>();
		nrBioAssemblies = -1;
	}

	/** String representation
	 *
	 */
	public String toString(){
		StringBuffer buf = new StringBuffer();

		try {

			@SuppressWarnings("rawtypes")
			Class c = Class.forName("org.biojava.bio.structure.PDBHeader");
			Method[] methods  = c.getMethods();

			for (int i = 0; i < methods.length; i++) {
				Method m = methods[i];

				String name = m.getName();

				if ( name.substring(0,3).equals("get")) {
					if (name.equals("getClass"))
						continue;
					Object o  = m.invoke(this, new Object[]{});
					if ( o != null){
						buf.append(name.substring(3,name.length()));
						buf.append(": " + o + " ");
						//if ( o instanceof Date) {
						//    buf.append(": " + FlatFileInstallation.dateFormat.format(o) + " ");
						//} else  {
						//    buf.append(": " + o + " ");
						//}
					}
				}
			}
		} catch (Exception e){
			e.printStackTrace();
		}

		return buf.toString();
	}

	/** Return a PDB representation of the PDB Header
	 *
	 * @return a PDB file style display
	 */
	public String toPDB(){
		StringBuffer buf = new StringBuffer();
		toPDB(buf);
		return buf.toString();
	}

	/** Appends a PDB representation of the PDB header to the provided StringBuffer
	 *
	 * @param buf
	 */
	public void toPDB(StringBuffer buf){
		//          1         2         3         4         5         6         7
		//01234567890123456789012345678901234567890123456789012345678901234567890123456789
		//HEADER    COMPLEX (SERINE PROTEASE/INHIBITORS)    06-FEB-98   1A4W
		//TITLE     CRYSTAL STRUCTURES OF THROMBIN WITH THIAZOLE-CONTAINING
		//TITLE    2 INHIBITORS: PROBES OF THE S1' BINDING SITE

		printHeader(buf);
		printTitle(buf);
		printExpdata(buf);
		printAuthors(buf);
		printResolution(buf);

	}

	private void printResolution(StringBuffer buf){

		if (getResolution() == DEFAULT_RESOLUTION){
			return;
		}

		DecimalFormat d2 = (DecimalFormat)NumberFormat.getInstance(java.util.Locale.UK);
		d2.setMaximumIntegerDigits(2);
		d2.setMinimumFractionDigits(2);
		d2.setMaximumFractionDigits(2);

		buf.append("REMARK   2 RESOLUTION. ");
		String x = d2.format(resolution);
		buf.append(x);
		buf.append(" ANGSTROMS.");
		fillLine(buf,34+x.length());

		buf.append(newline);
	}

	private void printExpdata(StringBuffer buf){
		String exp = getTechnique();
		if ( exp == null || exp.length() == 0)
			return;

		buf.append("EXPDTA    ");
		buf.append(exp);

		// fill up the white space to the right column
		int l =  exp.length()+ 10;
		fillLine(buf,l);

		buf.append(newline);

	}

	private void printAuthors(StringBuffer buf){
		String authors = getAuthors();
		if ( authors == null)
			return;
		if ( authors.equals("")){
			return;
		}

		printMultiLine(buf, "AUTHOR   ", authors,',');

	}

	private void printMultiLine(StringBuffer buf, String lineStart, String data, char breakChar){
		if ( lineStart.length() !=  9)
			System.err.println("lineStart != 9, there will be problems :" + lineStart);

		if ( data.length() < 58) {
			buf.append(lineStart);
			buf.append(" ");
			buf.append(data);
			buf.append(newline);
			return;
		}
		String thisLine = "";
		int count = 1;
		while (data.length() > 57) {
			// find first whitespace from left
			// there are 10 chars to the left, so the cutoff position is 56
			boolean charFound = false;
			for ( int i =57;i>-1;i--){
				char c = data.charAt(i);
				if (c == breakChar){                
					// found the whitespace

					thisLine = data.substring(0,i+1);

					// prevent endless loop
					if (i == 0 )
						i++;
					data = data.substring(i);
					charFound = true;
					//System.out.println(thisLine);
					//System.out.println(title);
					break;
				}
			}
			// for emergencies...  prevents an endless loop
			if ( ! charFound){
				thisLine = data.substring(0,58);
				data = data.substring(57);             
			}
			if ( ( breakChar == ',' ) && ( data.charAt(0)== ',')) {
				data =   data.substring(1);
			}

			//TODO: check structures that have more than 10  lines...
			// start printing..

			buf.append(lineStart);
			if ( count > 1) {
				buf.append(count);
				if ( breakChar != ' ' )
					buf.append(" ");
			}
			else
				buf.append(" ");
			buf.append(thisLine);

			// fill up the white space to the right column
			int l =  thisLine.length()+ 10;
			while (l < 67){
				l++;
				buf.append(" ");
			}

			buf.append(newline);
			count++;

		}

		// last line...
		if ( data.trim().length() > 0){
			buf.append(lineStart);
			buf.append(count);
			int filledLeft = 10;
			if ( breakChar != ' ' ) {
				buf.append(" ");
				filledLeft++;
			}
			buf.append(data);
			// fill up the white space to the right column
			int l =  data.length()+ filledLeft;
			fillLine(buf,l);
			buf.append(newline);
		}

	}

	private void fillLine(StringBuffer buf, int currentPos){
		int l = currentPos;
		while (l < 67){
			l++;
			buf.append(" ");
		}
	}

	private void printHeader(StringBuffer buf){

		String classification = getClassification();

		if (
				(classification ==null) ||
				(classification.length() == 0))
			return;

		// we can;t display this line since the classification is not there...

		buf.append("HEADER    ");
		buf.append(classification);
		buf.append(" ");

		// fill up the white space to the right column
		int l =  classification.length() + 10 ;
		while (l < 49){
			l++;
			buf.append(" ");
		}

		Date d = getDepDate();
		if ( d !=  null){
			// provide correct display of Dep date...
			buf.append(dateFormat.format(d));
		} else {
			buf.append("         ");
		}
		buf.append("   ");

		String id = getIdCode();
		if ( id != null){
			buf.append(getIdCode());
			buf.append(" ");
		}
		else
			buf.append("    ");
		buf.append(newline);


	}

	private void printTitle(StringBuffer buf) {
		//          1         2         3         4         5         6         7
		//01234567890123456789012345678901234567890123456789012345678901234567890123456789

		//HEADER    COMPLEX (SERINE PROTEASE/INHIBITORS)    06-FEB-98   1A4W
		//TITLE     CRYSTAL STRUCTURES OF THROMBIN WITH THIAZOLE-CONTAINING
		//TITLE    2 INHIBITORS: PROBES OF THE S1' BINDING SITE

		String title = getTitle();

		if ( (title == null) || ( title.trim().length() == 0) )
			return;

		printMultiLine(buf, "TITLE    ", title,' ');

	}

	/** Get the ID used by Hibernate.
	 *
	 * @return the ID used by Hibernate
	 * @see #setId(Long)
	 */
	public Long getId() {
		return id;
	}

	/** Set the ID used by Hibernate.
	 *
	 * @param id the id assigned by Hibernate
	 * @see #getId()
	 *
	 */

	@SuppressWarnings("unused")
	private void setId(Long id) {
		this.id = id;
	}

	/** Compare two PDBHeader objects
	 *
	 * @param other a PDBHeader object to compare this one to.
	 * @return true if they are equal or false if they are not.
	 */
	public boolean equals(PDBHeader other){
		try {

			@SuppressWarnings("rawtypes")
			Class c = Class.forName("org.biojava.bio.structure.PDBHeader");
			Method[] methods  = c.getMethods();

			for (int i = 0; i < methods.length; i++) {
				Method m = methods[i];
				String name = m.getName();

				if ( name.substring(0,3).equals("get")) {
					if (name.equals("getClass"))
						continue;
					Object a  = m.invoke(this,  new Object[]{});
					Object b  = m.invoke(other, new Object[]{});
					if ( a == null ){
						if ( b == null ){
							continue;
						} else {
							System.out.println(name + " a is null, where other is " + b);
							return false;
						}
					}
					if ( b == null) {
						System.out.println(name + " other is null, where a is " + a);
						return false;
					}
					if (! (a.equals(b))){
						System.out.println("mismatch with " + name + " >" + a + "< >" + b + "<");
						return false;
					}
				}
			}
		} catch (Exception e){
			e.printStackTrace();
			return false;
		}
		return true;
	}


	/** The PDB code for this protein structure.
	 *
	 * @return the PDB identifier
	 * @see #setIdCode(String)
	 */
	public String getIdCode() {
		return idCode;
	}
	/** The PDB code for this protein structure.
	 *
	 * @param idCode the PDB identifier
	 * @see #getIdCode()
	 *
	 */
	public void setIdCode(String idCode) {
		this.idCode = idCode;
	}


	public String getClassification() {
		return classification;
	}

	public void setClassification(String classification) {
		this.classification = classification;
	}

	public Date getDepDate() {
		return depDate;
	}

	public void setDepDate(Date depDate) {
		this.depDate = depDate;
	}

	public String getTechnique() {
		return technique;
	}

	public void setTechnique(String technique) {
		this.technique = technique;
	}

	public float getResolution() {
		return resolution;
	}

	public void setResolution(float resolution) {
		this.resolution = resolution;
	}

	public Date getModDate() {
		return modDate;
	}

	public void setModDate(Date modDate) {
		this.modDate = modDate;
	}

	public String getMethod() {
		return method;
	}
	public void setMethod(String method) {
		this.method = method;
	}
	public String getTitle() {
		return title;
	}
	public void setTitle(String title) {
		this.title = title;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}

	/** Returns the names of the authors as listed in the AUTHORS section of a PDB file.
	 * Not necessarily the same authors as listed in the AUTH section of the primary citation!
	 *
	 * @return Authors as a string
	 */
	public String getAuthors()
	{
		return authors;
	}

	public void setAuthors(String authors)
	{
		this.authors = authors;
	}

	public Map<Integer,List<ModelTransformationMatrix>> getBioUnitTranformationMap() {
		return tranformationMap ;
	}

	public void setBioUnitTranformationMap(Map<Integer,List<ModelTransformationMatrix>> tranformationMap) {
		this.tranformationMap = tranformationMap;
	}

	public void setNrBioAssemblies(int nrBioAssemblies) {
		this.nrBioAssemblies = nrBioAssemblies;
		
	}
	
	public int getNrBioAssemblies() {
		return nrBioAssemblies;
	}



}
