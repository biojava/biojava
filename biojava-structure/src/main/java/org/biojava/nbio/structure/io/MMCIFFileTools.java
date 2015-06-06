package org.biojava.nbio.structure.io;

import java.lang.reflect.Field;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MMCIFFileTools {

	private static final Logger logger = LoggerFactory.getLogger(MMCIFFileTools.class);
	
	private static final String newline = System.getProperty("line.separator");
	
	public static String toLoopMmCifHeaderString(String categoryName, String className) throws ClassNotFoundException {
		StringBuilder str = new StringBuilder();
		
		str.append("loop_"+newline);
		
		Class<?> c = Class.forName(className);
		
		for (Field f : c.getDeclaredFields()) {
			str.append(categoryName+"."+f.getName()+newline);
		}
		
		return str.toString();
	}
	
	public static String toSingleLineMmCifString(Object a) {
		
		StringBuilder str = new StringBuilder();
		
		Class<?> c = a.getClass();
		
		for (Field f : c.getDeclaredFields()) {
			f.setAccessible(true);

			try {
				Object obj = f.get(a);
				if (obj==null) {
					logger.info("Field {} is null, will not write it out",f.getName()); 
					continue;
				}
				String val = (String) obj;
				
				str.append(String.format("%-9s", addMmCifQuoting(val)));
								
				
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
		
		//return 
		//	String.format(
		//		"%-6s %-7s %-2s %-4s %-1s"+newline, 
		//		a.getGroup_PDB(), a.getId(), a.getType_symbol(), a.getLabel_atom_id(), a.getLabel_alt_id());
		
	}
	
	private static String addMmCifQuoting(String val) {
		String newval;
		
		if (val.contains("'")) {
			// double quoting for strings containing single quotes
			newval = "\""+val+"\"";
		} else if (val.contains(" ")) {
			// single quoting for stings containing spaces
			newval = "'"+val+"'";
		} else {
			if (val.contains(" ") && val.contains("'")) {
				// TODO deal with this case
				logger.warn("Value contains both spaces and single quotes, won't format it: {}",val);
			}
			newval = val;
		}
		return newval;
	}
	
//	private static void addMmCifAtomSiteHeader(StringBuilder str) {
//	str.append("loop_"+newline);
//	
//	String atomSiteCategory = "_atom_site";
//	
//	str.append(atomSiteCategory+".group_PDB"+newline);
//	str.append(atomSiteCategory+".id"+newline);
//
//	str.append(atomSiteCategory+".type_symbol"+newline);
//
//	str.append(atomSiteCategory+".label_atom_id+"+newline);
//	str.append(atomSiteCategory+".label_alt_id"+newline);
//	str.append(atomSiteCategory+".label_comp_id"+newline);
//
//	str.append(atomSiteCategory+".label_asym_id"+newline);
//	str.append(atomSiteCategory+".label_entity_id"+newline);
//	str.append(atomSiteCategory+".label_seq_id"+newline);
//
//	str.append(atomSiteCategory+".pdbx_PDB_ins_code"+newline);
//	str.append(atomSiteCategory+".Cartn_x"+newline);
//	str.append(atomSiteCategory+".Cartn_y"+newline);
//	str.append(atomSiteCategory+".Cartn_z"+newline);
//	str.append(atomSiteCategory+".occupancy"+newline);
//	str.append(atomSiteCategory+".B_iso_or_equiv"+newline);
//
//	str.append(atomSiteCategory+".pdbx_formal_charge"+newline);
//
//	str.append(atomSiteCategory+".auth_seq_id"+newline);
//	str.append(atomSiteCategory+".auth_asym_id"+newline);
//
//	str.append(atomSiteCategory+".pdbx_PDB_model_num"+newline);
//	
//}
}
