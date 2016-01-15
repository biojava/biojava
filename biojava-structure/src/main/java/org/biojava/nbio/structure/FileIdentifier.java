package org.biojava.nbio.structure;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.structure.align.util.AtomCache;

public class FileIdentifier implements StructureIdentifier {
	private static final Pattern PDB_REGEX = Pattern.compile("([0-9][a-z]{3})[._].*",Pattern.CASE_INSENSITIVE);

	private final File file;
	public FileIdentifier(File file) {
		this.file = file;
	}
	public FileIdentifier(String name) {
		this(new File(name));
	}
	
	public File getFile() {
		return file;
	}

	@Override
	public String getIdentifier() {
		return file.toString();
	}

	@Override
	public Structure loadStructure(AtomCache cache) throws StructureException,
			IOException {
		return null;
	}

	/**
	 * Represents the full substructure.
	 * 
	 * Attempts to guess the PDB ID from the filename, but may give up and set
	 * it to null.
	 */
	@Override
	public SubstructureIdentifier toCanonical() throws StructureException {
		return new SubstructureIdentifier(guessPDBID(file.getName()), new ArrayList<ResidueRange>());
	}

	private static String guessPDBID(String name) {
		Matcher match = PDB_REGEX.matcher(name);
		if(match.matches()) {
			return match.group(0);
		} else {
			// Give up if doesn't match
			return null;
		}
	}
	/**
	 * Returns the complete input structure
	 */
	@Override
	public Structure reduce(Structure input) throws StructureException {
		return input;
	}

}
