/*
 *                    BioJava development code
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
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.mmcif.model.DatabasePdbrevRecord;
import org.biojava.nbio.structure.quaternary.BioAssemblyInfo;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.*;


/** 
 * A class that contains PDB Header information.
 *
 * @author Andreas Prlic
 * @since 1.6
 *
 */
public class PDBHeader implements PDBRecord {

	private static final long serialVersionUID = -5834326174085429508L;

	private static final Logger logger = LoggerFactory.getLogger(PDBHeader.class);

	private String title;
	private String description;
	private String idCode;
	private String classification;

	private Date depDate;
	private Date relDate;
	private Date modDate;

	private Set<ExperimentalTechnique> techniques;
	private PDBCrystallographicInfo crystallographicInfo;

	private float resolution;
	private float rFree;
	private float rWork;

	private JournalArticle journalArticle;
	private String authors;

	public static final float DEFAULT_RESOLUTION = 99;
	public static final float DEFAULT_RFREE = 1; // worst possible rfree is the default


	private Long id;
	public static final String newline = System.getProperty("line.separator");

	private DateFormat dateFormat;

	private Map<Integer,BioAssemblyInfo> bioAssemblies ;

	List<DatabasePdbrevRecord> revisionRecords;

	public PDBHeader(){

		depDate = new Date(0);
		modDate = new Date(0);
		relDate = new Date(0);
		dateFormat = new SimpleDateFormat("dd-MMM-yy",Locale.US);
		
		resolution = DEFAULT_RESOLUTION;
		rFree = DEFAULT_RFREE;
		rWork = DEFAULT_RFREE;
		
		bioAssemblies = new LinkedHashMap<Integer, BioAssemblyInfo>();
		crystallographicInfo = new PDBCrystallographicInfo();

	}

	/** String representation
	 *
	 */
	@Override
	public String toString(){
		StringBuilder buf = new StringBuilder();

		try {


			Class<?> c = Class.forName(PDBHeader.class.getName());
			Method[] methods  = c.getMethods();

			for (Method m : methods) {
				String name = m.getName();

				if (name.substring(0, 3).equals("get")) {
					if (name.equals("getClass")) {
						continue;
					}
					Object o = m.invoke(this);
					if (o != null) {
						buf.append(name.substring(3, name.length()));
						buf.append(": ").append(o).append(" ");
					}
				}
			}
		} catch (ClassNotFoundException e) {
			logger.error("Exception caught while creating toString  ",e);
		} catch (InvocationTargetException e) {
			logger.error("Exception caught while creating toString ",e);
		} catch (IllegalAccessException e) {
			logger.error("Exception caught while creating toString ",e);
		}

		return buf.toString();
	}

	/** Return a PDB representation of the PDB Header
	 *
	 * @return a PDB file style display
	 */
	@Override
	public String toPDB(){
		StringBuffer buf = new StringBuffer();
		toPDB(buf);
		return buf.toString();
	}

	/** Appends a PDB representation of the PDB header to the provided StringBuffer
	 *
	 * @param buf
	 */
	@Override
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
		Set<ExperimentalTechnique> exp = getExperimentalTechniques();
		if ( exp == null )
			return;


		buf.append("EXPDTA    ");

		int length = 0;
		int i = 0;
		for (ExperimentalTechnique et:exp) {
			if (i>0) {
				buf.append("; ");
				length+=2;
			}
			buf.append(et.getName());
			length+=et.getName().length();
			i++;
		}

		// fill up the white space to the right column
		int l =  length + 10;
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
			logger.info("lineStart != 9, there will be problems :" + lineStart);

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
		if (!data.trim().isEmpty()){
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

		if (classification == null || classification.isEmpty()) return;

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

		if ( (title == null) || (title.trim().isEmpty()) )
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

			Class<?> c = Class.forName(PDBHeader.class.getName());
			Method[] methods  = c.getMethods();

			for (Method m : methods) {
				String name = m.getName();

				if (name.substring(0, 3).equals("get")) {
					if (name.equals("getClass")) {
						continue;
					}
					Object a = m.invoke(this);
					Object b = m.invoke(other);
					if (a == null) {
						if (b == null) {
							continue;
						} else {
							logger.warn(name + " a is null, where other is " + b);
							return false;
						}
					}
					if (b == null) {
						logger.warn(name + " other is null, where a is " + a);
						return false;
					}
					if (!(a.equals(b))) {
						logger.warn("mismatch with " + name + " >" + a + "< >" + b + "<");
						return false;
					}
				}
			}
		} catch (ClassNotFoundException e) {
			logger.error("Exception caught while comparing PDBHeader objects ",e);
			return false;
		} catch (InvocationTargetException e) {
			logger.error("Exception caught while comparing PDBHeader objects ",e);
			return false;
		} catch (IllegalAccessException e) {
			logger.error("Exception caught while comparing PDBHeader objects ",e);
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

	/**
	 * Return the deposition date of the structure in the PDB.
	 * 
	 * @return the deposition date
	 */
	public Date getDepDate() {
		return depDate;
	}

	/**
	 * The deposition date of the structure in the PDB
	 * 
	 * @param depDate the deposition date
	 */
	public void setDepDate(Date depDate) {
		this.depDate = depDate;
	}

	/**
	 * Return the Set of ExperimentalTechniques, usually the set is of size 1 except for hybrid
	 * experimental techniques when the Set will contain 2 or more values
	 * @return the Set of ExperimentalTechniques or null if not set
	 */
	public Set<ExperimentalTechnique> getExperimentalTechniques() {
		return techniques;
	}

	/**
	 * Adds the experimental technique to the set of experimental techniques of this header.
	 * Note that if input is not a recognised technique string then no errors will be produced but
	 * false will be returned
	 * @param techniqueStr
	 * @return true if the input corresponds to a recognised technique string (see {@link ExperimentalTechnique})
	 * and it was not already present in the current set of ExperimentalTechniques
	 */
	public boolean setExperimentalTechnique(String techniqueStr) {

		ExperimentalTechnique et = ExperimentalTechnique.getByName(techniqueStr);

		if (et==null) return false;

		if (techniques==null) {
			techniques = EnumSet.of(et);
			return true;
		} else {
			return techniques.add(et);
		}

	}

	public PDBCrystallographicInfo getCrystallographicInfo() {
		return crystallographicInfo;
	}

	public void setCrystallographicInfo(PDBCrystallographicInfo crystallographicInfo) {
		this.crystallographicInfo = crystallographicInfo;
	}

	public float getResolution() {
		return resolution;
	}

	public void setResolution(float resolution) {
		this.resolution = resolution;
	}

	public float getRfree() {
		return rFree;
	}

	public void setRfree(float rFree) {
		this.rFree = rFree;
	}

	/**
	 * Return the latest modification date of the structure.
	 * 
	 * @return the latest modification date
	 */
	public Date getModDate() {
		return modDate;
	}

	/**
	 * The latest modification date of the structure.
	 * 
	 * @param modDate the latest modification date
	 */
	public void setModDate(Date modDate) {
		this.modDate = modDate;
	}
	
	/**
	 * Return the release date of the structure in the PDB.
	 * 
	 * @return the release date
	 */
	public Date getRelDate() {
		return relDate;
	}

	/**
	 * 
	 * The release date of the structure in the PDB.
	 * 
	 * @param relDate the release date
	 */
	public void setRelDate(Date relDate) {
		this.relDate = relDate;
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

	/**
	 * Return the names of the authors as listed in the AUTHORS section of a PDB file.
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

	/**
	 * Return whether or not the entry has an associated journal article
	 * or publication. The JRNL section is not mandatory and thus may not be
	 * present.
	 * @return flag if a JournalArticle could be found.
	 */
	public boolean hasJournalArticle() {
		return this.journalArticle != null;
	}

	/**
	 * Get the associated publication as defined by the JRNL records in a PDB
	 * file.
	 * @return a JournalArticle
	 */
	public JournalArticle getJournalArticle() {
		return this.journalArticle;
	}

	/**
	 * Set the associated publication as defined by the JRNL records in a PDB
	 * file.
	 * @param journalArticle the article
	 */
	public void setJournalArticle(JournalArticle journalArticle) {
		this.journalArticle = journalArticle;
	}

	/**
	 * Return the map of biological assemblies. The keys are the
	 * biological assembly identifiers (starting at 1). Non-numerical identifiers
	 * such as PAU or XAU are not supported.
	 * @return
	 */
	public Map<Integer,BioAssemblyInfo> getBioAssemblies() {
		return bioAssemblies ;
	}

	public void setBioAssemblies(Map<Integer,BioAssemblyInfo> bioAssemblies) {
		this.bioAssemblies = bioAssemblies;
	}

	/**
	 * Get the number of biological assemblies available in the PDB header
	 * @return
	 */
	public int getNrBioAssemblies() {
		return this.bioAssemblies.size();
	}

	public List<DatabasePdbrevRecord> getRevisionRecords() {
		return revisionRecords;
	}

	public void setRevisionRecords(List<DatabasePdbrevRecord> revisionRecords) {
		this.revisionRecords = revisionRecords;
	}

	/**
	 * @return the R-work for this structure.
	 */
	public float getRwork() {
		return rWork;
	}

	/**
	 * @param rWork  the R-work for this structure.
	 */
	public void setRwork(float rWork) {
		this.rWork = rWork;
	}
}
