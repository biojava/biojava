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
package org.biojava.nbio.structure.io.mmcif;


import java.lang.reflect.Field;
import java.util.*;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.FileConvert;
import org.biojava.nbio.structure.io.mmcif.model.AbstractBean;
import org.biojava.nbio.structure.io.mmcif.model.AtomSite;
import org.biojava.nbio.structure.io.mmcif.model.CIFLabel;
import org.biojava.nbio.structure.io.mmcif.model.Cell;
import org.biojava.nbio.structure.io.mmcif.model.IgnoreField;
import org.biojava.nbio.structure.io.mmcif.model.Symmetry;
import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.SpaceGroup;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Some tools for mmCIF file writing.
 *
 * See http://www.iucr.org/__data/assets/pdf_file/0019/22618/cifguide.pdf
 *
 * CIF categories are represented as a simple bean, typically extending {@link AbstractBean}.
 * By default, all fields from the bean are taken as the CIF labels. Fields
 * may be omitted by annotating them as {@link IgnoreField @IgnoreField}.
 * The CIF label for a field may be changed (for instance, for fields that
 * are not valid Java identifiers) by defining a function
 * <tt>static Map<String,String> getCIFLabelMap()</tt>
 * mapping from the field's name to the correct label.
 * 
 * @author Jose Duarte
 * @author Spencer Bliven
 */
public class MMCIFFileTools {

	private static final Logger logger = LoggerFactory.getLogger(MMCIFFileTools.class);

	private static final String newline = System.getProperty("line.separator");

	/**
	 * The character to be printed out in cases where a value is not assigned in mmCIF files
	 */
	public static final String MMCIF_MISSING_VALUE = "?";

	/**
	 * The character to be printed out as a default value in mmCIF files, e.g. for the default alt_locs
	 */
	public static final String MMCIF_DEFAULT_VALUE = ".";


	/**
	 * Produces a mmCIF loop header string for the given categoryName and className.
	 * className must be one of the beans in the {@link org.biojava.nbio.structure.io.mmcif.model} package
	 * @param categoryName
	 * @param className
	 * @return
	 * @throws ClassNotFoundException if the given className can not be found
	 */
	public static String toLoopMmCifHeaderString(String categoryName, String className) throws ClassNotFoundException {
		StringBuilder str = new StringBuilder();

		str.append(SimpleMMcifParser.LOOP_START+newline);

		Class<?> c = Class.forName(className);

		for (Field f : getFields(c)) {
			str.append(categoryName+"."+f.getName()+newline);
		}

		return str.toString();
	}

	/**
	 * Converts a mmCIF bean (see {@link org.biojava.nbio.structure.io.mmcif.model} to
	 * a String representing it in mmCIF (single-record) format.
	 * @param categoryName
	 * @param o
	 * @return
	 */
	public static String toMMCIF(String categoryName, Object o) {

		StringBuilder sb = new StringBuilder();

		Class<?> c = o.getClass();


		Field[] fields = getFields(c);
		String[] names = getFieldNames(fields);

		int maxFieldNameLength = getMaxStringLength(names);

		for (int i=0;i<fields.length;i++) {
			Field f = fields[i];
			String name = names[i];

			sb.append(categoryName).append(".").append(name);

			int spacing = maxFieldNameLength - name.length() + 3;

			try {
				Object obj = f.get(o);
				String val;
				if (obj==null) {
					logger.debug("Field {} is null, will write it out as {}",name,MMCIF_MISSING_VALUE);
					val = MMCIF_MISSING_VALUE;
				} else {
					val = (String) obj;
				}
				for (int j=0;j<spacing;j++) sb.append(' ');
				sb.append(addMmCifQuoting(val));
				sb.append(newline);

			} catch (IllegalAccessException e) {
				logger.warn("Field {} is inaccessible", name);
				continue;
			} catch (ClassCastException e) {
				logger.warn("Could not cast value to String for field {}",name);
				continue;
			}

		}

		sb.append(SimpleMMcifParser.COMMENT_CHAR+newline);

		return sb.toString();
	}

	/**
	 * Gets all fields for a particular class, filtering fields annotated
	 * with {@link IgnoreField @IgnoreField}.
	 * 
	 * As a side effect, calls {@link Field#setAccessible(boolean) setAccessible(true)}
	 * on all fields.
	 * @param c
	 * @return
	 */
	public static Field[] getFields(Class<?> c) {
		Field[] allFields = c.getDeclaredFields();
		Field[] fields = new Field[allFields.length];
		int n = 0;
		for(Field f : allFields) {
			f.setAccessible(true);
			IgnoreField anno = f.getAnnotation(IgnoreField.class);
			if(anno == null) {
				fields[n] = f;
				n++;
			}
		}
		return Arrays.copyOf(fields, n);
	}

	/**
	 * Gets the mmCIF record name for each field. This is generally just
	 * the name of the field or the value specified by the {@link CIFLabel @CIFLabel} annotation.
	 * 
	 * As a side effect, calls {@link Field#setAccessible(boolean) setAccessible(true)}
	 * on all fields.
	 * @param fields
	 * @return
	 */
	public static String[] getFieldNames(Field[] fields) {
		String[] names = new String[fields.length];
		for(int i=0;i<fields.length;i++) {
			Field f = fields[i];
			f.setAccessible(true);
			String rawName = fields[i].getName();
			CIFLabel cifLabel = f.getAnnotation(CIFLabel.class);
			if(cifLabel != null) {
				names[i] = cifLabel.label();
			} else {
				names[i] = rawName;
			}
		}
		return names;
	}

	/**
	 * Converts a list of mmCIF beans (see {@link org.biojava.nbio.structure.io.mmcif.model} to
	 * a String representing them in mmCIF loop format with one record per line.
	 * @param list
	 * @return
	 */
	public static <T> String toMMCIF(List<T> list, Class<T> klass) {
		if (list.isEmpty()) throw new IllegalArgumentException("List of beans is empty!");

		Field[] fields = getFields(klass);
		int[] sizes = getFieldSizes(list,fields);

		StringBuilder sb = new StringBuilder();

		for (T o:list) {
			sb.append(toSingleLoopLineMmCifString(o, fields, sizes));
		}

		sb.append(SimpleMMcifParser.COMMENT_CHAR+newline);

		return sb.toString();
	}

	/**
	 * Given a mmCIF bean produces a String representing it in mmCIF loop format as a single record line
	 * @param record
	 * @param fields Set of fields for the record. If null, will be calculated from the class of the record
	 * @param sizes the size of each of the fields
	 * @return
	 */
	private static String toSingleLoopLineMmCifString(Object record, Field[] fields, int[] sizes) {

		StringBuilder str = new StringBuilder();

		Class<?> c = record.getClass();

		if(fields == null)
			fields = getFields(c);
		
		if (sizes.length!=fields.length)
			throw new IllegalArgumentException("The given sizes of fields differ from the number of declared fields");

		int i = -1;
		for (Field f : fields) {
			i++;
			f.setAccessible(true);

			try {
				Object obj = f.get(record);
				String val;
				if (obj==null) {
					logger.debug("Field {} is null, will write it out as {}",f.getName(),MMCIF_MISSING_VALUE);
					val = MMCIF_MISSING_VALUE;
				} else {
					val = (String) obj;
				}

				str.append(String.format("%-"+sizes[i]+"s ", addMmCifQuoting(val)));


			} catch (IllegalAccessException e) {
				logger.warn("Field {} is inaccessible", f.getName());
				continue;
			} catch (ClassCastException e) {
				logger.warn("Could not cast value to String for field {}",f.getName());
				continue;
			}
		}

		str.append(newline);

		return str.toString();

	}

	/**
	 * Adds quoting to a String according to the STAR format (mmCIF) rules
	 * @param val
	 * @return
	 */
	private static String addMmCifQuoting(String val) {
		String newval;

		if (val.contains("'")) {
			// double quoting for strings containing single quotes (not strictly necessary but it's what the PDB usually does)
			newval = "\""+val+"\"";
		} else if (val.contains(" ")) {
			// single quoting for stings containing spaces
			newval = "'"+val+"'";
		} else {
			if (val.contains(" ") && val.contains("'")) {
				// TODO deal with this case
				logger.warn("Value contains both spaces and single quotes, won't format it: {}. CIF ouptut will likely be invalid.",val);
			}
			newval = val;
		}
		// TODO deal with all the other cases: e.g. multi-line quoting with ;;

		return newval;
	}

	/**
	 * Converts a SpaceGroup object to a {@link Symmetry} object.
	 * @param sg
	 * @return
	 */
	public static Symmetry convertSpaceGroupToSymmetry(SpaceGroup sg) {
		Symmetry sym = new Symmetry();
		sym.setSpace_group_name_H_M(sg.getShortSymbol());
		// TODO do we need to fill any of the other values?
		return sym;
	}

	/**
	 * Converts a CrystalCell object to a {@link Cell} object.
	 * @param c
	 * @return
	 */
	public static Cell convertCrystalCellToCell(CrystalCell c) {
		Cell cell = new Cell();
		cell.setLength_a(String.format("%.3f",c.getA()));
		cell.setLength_b(String.format("%.3f",c.getB()));
		cell.setLength_c(String.format("%.3f",c.getC()));
		cell.setAngle_alpha(String.format("%.3f",c.getAlpha()));
		cell.setAngle_beta(String.format("%.3f",c.getBeta()));
		cell.setAngle_gamma(String.format("%.3f",c.getGamma()));

		return cell;
	}

	/**
	 * Converts an Atom object to an {@link AtomSite} object.
	 * @param a
	 * @param model the model number for the output AtomSites
	 * @param chainName the chain identifier (author id) for the output AtomSites
	 * @param chainId the internal chain identifier (asym id) for the output AtomSites
	 * @return
	 */
	public static AtomSite convertAtomToAtomSite(Atom a, int model, String chainName, String chainId) {
		return convertAtomToAtomSite(a, model, chainName, chainId, a.getPDBserial());
	}

	/**
	 * Converts an Atom object to an {@link AtomSite} object.
	 * @param a the atom
	 * @param model the model number for the output AtomSites
	 * @param chainName the chain identifier (author id) for the output AtomSites
	 * @param chainId the internal chain identifier (asym id) for the output AtomSites
	 * @param atomId the atom id to be written to AtomSite
	 * @return
	 */
	public static AtomSite convertAtomToAtomSite(Atom a, int model, String chainName, String chainId, int atomId) {

		/*
		ATOM 7    C CD  . GLU A 1 24  ? -10.109 15.374 38.853 1.00 50.05 ? ? ? ? ? ? 24  GLU A CD  1
		ATOM 8    O OE1 . GLU A 1 24  ? -9.659  14.764 37.849 1.00 49.80 ? ? ? ? ? ? 24  GLU A OE1 1
		ATOM 9    O OE2 . GLU A 1 24  ? -11.259 15.171 39.310 1.00 50.51 ? ? ? ? ? ? 24  GLU A OE2 1
		ATOM 10   N N   . LEU A 1 25  ? -5.907  18.743 37.412 1.00 41.55 ? ? ? ? ? ? 25  LEU A N   1
		ATOM 11   C CA  . LEU A 1 25  ? -5.168  19.939 37.026 1.00 37.55 ? ? ? ? ? ? 25  LEU A CA  1
		*/

		Group g = a.getGroup();

		String record ;
		if ( g.getType().equals(GroupType.HETATM) ) {
			record = "HETATM";
		} else {
			record = "ATOM";
		}

		String entityId = "0";
		String labelSeqId = Integer.toString(g.getResidueNumber().getSeqNum());
		if (g.getChain()!=null && g.getChain().getEntityInfo()!=null) {
			entityId = Integer.toString(g.getChain().getEntityInfo().getMolId());
			labelSeqId = Integer.toString(g.getChain().getEntityInfo().getAlignedResIndex(g, g.getChain()));
		}

		Character  altLoc = a.getAltLoc()           ;
		String altLocStr;
		if (altLoc==null || altLoc == ' ') {
			altLocStr = MMCIF_DEFAULT_VALUE;
		} else {
			altLocStr = altLoc.toString();
		}

		Element e = a.getElement();
		String eString = e.toString().toUpperCase();
		if ( e.equals(Element.R)) {
			eString = "X";
		}

		String insCode = MMCIF_MISSING_VALUE;
		if (g.getResidueNumber().getInsCode()!=null ) {
			insCode = Character.toString(g.getResidueNumber().getInsCode());
		}

		AtomSite atomSite = new AtomSite();
		atomSite.setGroup_PDB(record);
		atomSite.setId(Integer.toString(atomId));
		atomSite.setType_symbol(eString);
		atomSite.setLabel_atom_id(a.getName());
		atomSite.setLabel_alt_id(altLocStr);
		atomSite.setLabel_comp_id(g.getPDBName());
		atomSite.setLabel_asym_id(chainId);
		atomSite.setLabel_entity_id(entityId);
		atomSite.setLabel_seq_id(labelSeqId);
		atomSite.setPdbx_PDB_ins_code(insCode);
		atomSite.setCartn_x(FileConvert.d3.format(a.getX()));
		atomSite.setCartn_y(FileConvert.d3.format(a.getY()));
		atomSite.setCartn_z(FileConvert.d3.format(a.getZ()));
		atomSite.setOccupancy(FileConvert.d2.format(a.getOccupancy()));
		atomSite.setB_iso_or_equiv(FileConvert.d2.format(a.getTempFactor()));
		atomSite.setAuth_seq_id(Integer.toString(g.getResidueNumber().getSeqNum()));
		atomSite.setAuth_comp_id(g.getPDBName());
		atomSite.setAuth_asym_id(chainName);
		atomSite.setAuth_atom_id(a.getName());
		atomSite.setPdbx_PDB_model_num(Integer.toString(model));

		return atomSite;
	}

	/**
	 * Converts a Group into a List of {@link AtomSite} objects.
	 * Atoms in other altloc groups (different from the main group) are also included, removing possible duplicates
	 * via using the atom identifier to assess uniqueness.
	 * @param g the group
	 * @param model the model number for the output AtomSites
	 * @param chainName the chain identifier (author id) for the output AtomSites
	 * @param chainId the internal chain identifier (asym id) for the output AtomSites
	 * @return
	 */
	public static List<AtomSite> convertGroupToAtomSites(Group g, int model, String chainName, String chainId) {

		// The alt locs can have duplicates, since at parsing time we make sure that all alt loc groups have
		// all atoms (see StructureTools#cleanUpAltLocs)
		// Thus we have to remove duplicates here by using the atom id
		// See issue https://github.com/biojava/biojava/issues/778 and TestAltLocs.testMmcifWritingAllAltlocs/testMmcifWritingPartialAltlocs
		Map<Integer, AtomSite> uniqueAtomSites = new LinkedHashMap<>();

		int groupsize  = g.size();

		for ( int atompos = 0 ; atompos < groupsize; atompos++) {
			Atom a = g.getAtom(atompos);
			if ( a == null)
				continue ;

			uniqueAtomSites.put(a.getPDBserial(), convertAtomToAtomSite(a, model, chainName, chainId));
		}

		if ( g.hasAltLoc()){
			for (Group alt : g.getAltLocs() ) {
				for (AtomSite atomSite : convertGroupToAtomSites(alt, model, chainName, chainId)) {
					uniqueAtomSites.put(Integer.parseInt(atomSite.getId()), atomSite);
				}
			}
		}
		return new ArrayList<>(uniqueAtomSites.values());
	}

	/**
	 * Converts a Chain into a List of {@link AtomSite} objects
	 * @param c the chain
	 * @param model the model number for the output AtomSites
	 * @param chainName the chain identifier (author id) for the output AtomSites
	 * @param chainId the internal chain identifier (asym id) for the output AtomSites
	 * @return
	 */
	public static List<AtomSite> convertChainToAtomSites(Chain c, int model, String chainName, String chainId) {

		List<AtomSite> list = new ArrayList<>();

		if (c.getEntityInfo()==null) {
			logger.warn("No Compound (entity) found for chain {}: entity_id will be set to 0, label_seq_id will be the same as auth_seq_id", c.getName());
		}

		for ( int h=0; h<c.getAtomLength();h++){

			Group g= c.getAtomGroup(h);

			list.addAll(convertGroupToAtomSites(g, model, chainName, chainId));

		}

		return list;
	}

	/**
	 * Converts a Structure into a List of {@link AtomSite} objects
	 * @param s
	 * @return
	 */
	public static List<AtomSite> convertStructureToAtomSites(Structure s) {
		List<AtomSite> list = new ArrayList<AtomSite>();

		for (int m=0;m<s.nrModels();m++) {
			for (Chain c:s.getChains(m)) {
				list.addAll(convertChainToAtomSites(c, m+1, c.getName(), c.getId()));
			}
		}
		return list;
	}

	/**
	 * Finds the max length of each of the String values contained in each of the fields of the given list of beans.
	 * Useful for producing mmCIF loop data that is aligned for all columns.
	 * @param list list of objects. All objects should have the same class.
	 * @param fields Set of fields for the record. If null, will be calculated from the class of the first record
	 * @return
	 * @see #toMMCIF(List, Class)
	 */
	private static <T> int[] getFieldSizes(List<T> list, Field[] fields) {

		if (list.isEmpty()) throw new IllegalArgumentException("List of beans is empty!");

		if(fields == null)
			fields = getFields(list.get(0).getClass());

		int[] sizes = new int [fields.length];


		for (T a:list) {
			int i = -1;
			for (Field f : fields) {
				i++;

				f.setAccessible(true);

				try {
					Object obj = f.get(a);
					int length;
					if (obj==null) {
						length = MMCIF_MISSING_VALUE.length();
					} else {
						String val = (String) obj;
						length = addMmCifQuoting(val).length();
					}

					if (length>sizes[i]) sizes[i] = length;

				} catch (IllegalAccessException e) {
					logger.warn("Field {} is inaccessible", f.getName());
					continue;
				} catch (ClassCastException e) {
					logger.warn("Could not cast value to String for field {}",f.getName());
					continue;
				}
			}
		}
		return sizes;
	}

	/**
	 * Finds the max length of a list of strings
	 * Useful for producing mmCIF single-record data that is aligned for all values.
	 * @param names
	 * @return
	 * @see #toMMCIF(String, Object)
	 */
	private static int getMaxStringLength(String[] names) {
		int size = 0;
		for(String s : names) {
			if(s.length()>size) {
				size = s.length();
			}
		}
		return size;
	}
}
