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
/*
 * Created on 2005-08-01
 */

package org.biojava.bio.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

/**
 * <p>
 * This object is able to read a substitution matrix file and constructs a short
 * matrix in memory. Every single element of the matrix can be accessed by the
 * method <code>getValueAt</code> with the parameters being two BioJava symbols.
 * This is why it is not necessary to access the matrix directly. If there is no
 * value for the two specified <code>Symbol</code>s an <code>Exception</code> is
 * thrown.
 * </p>
 * <p>
 * Substitution matrix files, are available at <a
 * href="ftp://ftp.ncbi.nlm.nih.gov/blast/matrices/"> the NCBI FTP
 * directory</a>.
 * </p>
 * 
 * @author Andreas Dr&auml;ger <andreas.draeger@uni-tuebingen.de>
 */
public class SubstitutionMatrix {
	
	/**
	 * 
	 */
	private Map<Symbol, Integer> rowSymbols, colSymbols;

	/**
	 * 
	 */
	private short[][] matrix;

	/**
	 * Minimal and maximal entry in this matrix
	 */
	private short min, max;

	/**
	 * The alphabet used by this matrix.
	 */
	private FiniteAlphabet alphabet;

	/**
	 * Name and description of this matrix.
	 */
	private String description, name;

	/**
	 * Just the new line symbol of the system.
	 */
	private static final String newLine = System.getProperty("line.separator");

	/**
	 * This constructs a <code>SubstitutionMatrix</code> object that contains
	 * two <code>Map</code> data structures having BioJava symbols as keys and
	 * the value being the index of the matrix containing the substitution
	 * score.
	 * 
	 * @param alpha
	 *            the alphabet of the matrix (e.g., DNA, RNA or PROTEIN, or
	 *            PROTEIN-TERM)
	 * @param matrixFile
	 *            the file containing the substitution matrix. Lines starting
	 *            with '<code>#</code>' are comments. The line starting with a
	 *            white space, is the table head. Every line has to start with
	 *            the one letter representation of the Symbol and then the
	 *            values for the exchange.
	 * @throws IOException
	 * @throws BioException
	 * @throws NumberFormatException
	 */
	public SubstitutionMatrix(FiniteAlphabet alpha, File matrixFile)
			throws BioException, NumberFormatException, IOException {
		this.alphabet = alpha;
		this.description = "";
		this.name = matrixFile.getName();
		this.rowSymbols = new HashMap<Symbol, Integer>();
		this.colSymbols = new HashMap<Symbol, Integer>();
		this.matrix = this.parseMatrix(matrixFile);
	}

	/**
	 * With this constructor it is possible to construct a SubstitutionMatrix
	 * object from a substitution matrix file. The given String contains a
	 * number of lines separated by
	 * <code>System.getProperty("line.separator")</code>. Everything else is the
	 * same than for the constructor above.
	 * 
	 * @param alpha
	 *            The <code>FiniteAlphabet</code> to use
	 * @param matrixString
	 * @param name
	 *            of the matrix.
	 * @throws BioException
	 * @throws IOException
	 * @throws NumberFormatException
	 */
	public SubstitutionMatrix(FiniteAlphabet alpha, String matrixString,
			String name) throws BioException, NumberFormatException,
			IOException {
		this.alphabet = alpha;
		this.description = "";
		this.name = name;
		this.rowSymbols = new HashMap<Symbol, Integer>();
		this.colSymbols = new HashMap<Symbol, Integer>();
		this.matrix = this.parseMatrix(matrixString);
		// this.printMatrix();
	}

	/**
	 * Constructs a SubstitutionMatrix with every Match and every Replace having
	 * the same expenses given by the parameters. Ambiguous symbols are not
	 * considered because there might be to many of them (for proteins).
	 * 
	 * @param alpha
	 * @param match
	 * @param replace
	 */
	public SubstitutionMatrix(FiniteAlphabet alpha, short match, short replace) {
		int i = 0, j = 0;

		this.alphabet = alpha;
		this.description = "Identity matrix. All replaces and all matches are treated equally.";
		this.name = "IDENTITY_" + match + "_" + replace;
		this.rowSymbols = new HashMap<Symbol, Integer>();
		this.colSymbols = new HashMap<Symbol, Integer>();
		this.matrix = new short[alpha.size()][alpha.size()];

		Symbol[] sym = new Symbol[alpha.size()];
		Iterator<Symbol> iter = alpha.iterator();

		for (i = 0; iter.hasNext(); i++) {
			sym[i] = iter.next();
			rowSymbols.put(sym[i], new Integer(i));
			colSymbols.put(sym[i], new Integer(i));
		}

		for (i = 0; i < alphabet.size(); i++)
			for (j = 0; j < alphabet.size(); j++)
				if (sym[i].getMatches().contains(sym[j]))
					matrix[i][j] = match;
				else
					matrix[i][j] = replace;

		// this.printMatrix();
	}

	/**
	 * This constructor can be used to guess the alphabet of this substitution
	 * matrix. However, it is recommended to apply another constructor if the
	 * alphabet is known.
	 * 
	 * @param file
	 *            A file containing a substitution matrix.
	 * @throws NumberFormatException
	 * @throws NoSuchElementException
	 * @throws BioException
	 * @throws IOException
	 */
	public SubstitutionMatrix(File file) throws NumberFormatException,
			NoSuchElementException, BioException, IOException {
		this(guessAlphabet(file), file);
	}

	/**
	 * This constructor can be used to guess the alphabet of this substitution
	 * matrix. However, it is recommended to apply another constructor if the
	 * alphabet is known.
	 * 
	 * @param reader
	 * @throws NumberFormatException
	 * @throws BioException
	 * @throws IOException
	 */
	public static SubstitutionMatrix getSubstitutionMatrix(BufferedReader reader)
			throws NumberFormatException, BioException, IOException {
		StringBuffer stringMatrix = new StringBuffer("");
		while (reader.ready()) {
			stringMatrix.append(reader.readLine());
			stringMatrix.append(newLine);
		}
		reader.close();
		String mat = stringMatrix.toString();
		FiniteAlphabet alpha = guessAlphabet(new BufferedReader(
				new StringReader(mat)));
		SubstitutionMatrix matrix = new SubstitutionMatrix(alpha, mat,
				"unknown");
		return matrix;
	}

	/**
	 * This method tries to identify the alphabet within a matrix file. This is
	 * necessary in cases where we do not know if this is a matrix for DNA, RNA
	 * or PROTEIN/PROTEIN-TERM.
	 * 
	 * @param file
	 * @return
	 * @throws IOException
	 * @throws BioException
	 * @throws NoSuchElementException
	 * @throws BioException
	 */
	private static FiniteAlphabet guessAlphabet(File file) throws IOException,
			NoSuchElementException, BioException {
		String fileName = file.getName().toLowerCase();
		if (fileName.contains("pam") || fileName.contains("blosum"))
			return (FiniteAlphabet) AlphabetManager
					.alphabetForName("PROTEIN-TERM");
		return guessAlphabet(new BufferedReader(new FileReader(file)));
	}

	/**
	 * This method guesses the alphabet of the given substituttion matrix which
	 * is required for the parser.
	 * 
	 * @param reader
	 * @return
	 * @throws IOException
	 * @throws BioException
	 */
	private static FiniteAlphabet guessAlphabet(BufferedReader reader)
			throws IOException, BioException {
		String line, trim;
		FiniteAlphabet alphabet = null;
		while (reader.ready()) {
			line = reader.readLine();
			if (line == null)
				break;
			trim = line.trim();
			if (trim.charAt(0) == '#')
				continue;
			else if ((line.charAt(0) == ' ') || (line.charAt(0) == '\t')) {
				String alphabets[] = new String[] { "DNA", "RNA", "PROTEIN",
						"PROTEIN-TERM" };
				SymbolTokenization symtok;
				for (int i = 0; i < alphabets.length; i++) {
					alphabet = (FiniteAlphabet) AlphabetManager
							.alphabetForName(alphabets[i]);
					symtok = alphabet.getTokenization("token");
					StringTokenizer st = new StringTokenizer(trim);
					boolean noError = true;
					for (int j = 0; st.hasMoreElements(); j++)
						try {
							symtok.parseToken(st.nextElement().toString());
						} catch (IllegalSymbolException exc) {
							noError = false;
							break;
						}
					if (noError)
						return alphabet;
				}
			}
		}
		throw new BioException(
				"Unknow alphabet used in this substitution matrix");
	}

	/**
	 * Reads a String representing the contents of a substitution matrix file.
	 * 
	 * @param matrixObj
	 * @return matrix
	 * @throws BioException
	 * @throws IOException
	 * @throws NumberFormatException
	 */
	private short[][] parseMatrix(Object matrixObj) throws BioException,
			NumberFormatException, IOException {
		int j = 0, rows = 0, cols = 0;
		SymbolTokenization symtok = alphabet.getTokenization("token");
		StringTokenizer st;
		String line, trim;

		this.min = Short.MAX_VALUE;
		this.max = Short.MIN_VALUE;
		/*
		 * First: count how many elements are in the matrix fill lines and rows
		 */
		Reader reader;
		if (matrixObj instanceof File)
			reader = new FileReader((File) matrixObj);
		else if (matrixObj instanceof String)
			reader = new StringReader(matrixObj.toString());
		else
			return null;
		BufferedReader br = new BufferedReader(reader);

		while (br.ready()) {
			line = br.readLine();
			if (line == null)
				break;
			trim = line.trim();
			if (trim.length() == 0)
				continue;
			if (trim.charAt(0) == '#') {
				description += line.substring(1);
				continue;
			} else if (!trim.startsWith(newLine)) {
				if ((line.charAt(0) == ' ') || (line.charAt(0) == '\t')) {
					st = new StringTokenizer(trim);
					for (j = 0; st.hasMoreElements(); j++) {
						colSymbols.put(symtok.parseToken(st.nextElement()
								.toString()), Integer.valueOf(j));
					}
					cols = j;
				} else {
					// the matrix.
					st = new StringTokenizer(line);
					if (st.hasMoreElements())
						rowSymbols.put(symtok.parseToken(st.nextElement()
								.toString()), Integer.valueOf(rows++));
				}
			}
		}
		br.close();

		short[][] matrix = new short[rows][cols];

		rows = 0;
		if (matrixObj instanceof File)
			reader = new FileReader((File) matrixObj);
		else if (matrixObj instanceof String)
			reader = new StringReader(matrixObj.toString());
		else
			return null;
		br = new BufferedReader(reader);

		/*
		 * Second reading. Fill the matrix.
		 */
		while (br.ready()) {
			line = br.readLine();
			if (line == null)
				break;
			trim = line.trim();
			if (trim.length() == 0 || trim.charAt(0) == '#')
				continue;
			else if ((line.charAt(0) == ' ') || (line.charAt(0) == '\t'))
				continue;
			else if (!trim.startsWith(newLine)) { // lines:
				st = new StringTokenizer(trim);
				if (st.hasMoreElements())
					st.nextElement(); // throw away Symbol at
				// beginning.
				for (j = 0; st.hasMoreElements(); j++) {// cols:
					matrix[rows][j] = (short) Math.round(Double.parseDouble(st
							.nextElement().toString()));
					if (matrix[rows][j] > max)
						max = matrix[rows][j]; // maximum.
					if (matrix[rows][j] < min)
						min = matrix[rows][j]; // minimum.
				}
				rows++;
			}
		}
		br.close();

		return matrix;
	}

	/**
	 * There are some substitution matrices containing more columns than lines.
	 * This has to do with the ambiguous symbols. Lines are always good, columns
	 * might not contain the whole information. The matrix is supposed to be
	 * symmetric anyway, so you can always set the ambiguous symbol to be the
	 * first argument.
	 * 
	 * @param row
	 *            Symbol of the line
	 * @param col
	 *            Symbol of the column
	 * @return expenses for the exchange of symbol row and symbol column.
	 * @throws BioException
	 */
	public short getValueAt(Symbol row, Symbol col) throws BioException {
		if ((!rowSymbols.containsKey(row)) || (!colSymbols.containsKey(col))) {
			System.err.printf("SubstitutionMatrix: No entry for the symbols %s and %s\n",
					row.getName(), col.getName());

			// treat the two records as X:
			return 0;
		}
		return matrix[rowSymbols.get(row).intValue()][colSymbols.get(col)
				.intValue()];
	}

	/**
	 * This gives you the description of this matrix if there is one. Normally
	 * substitution matrix files like BLOSUM contain some lines of description.
	 * 
	 * @return the comment of the matrix
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * Every substitution matrix has a name like "BLOSUM30" or "PAM160". This
	 * will be returned by this method.
	 * 
	 * @return the name of the matrix.
	 */
	public String getName() {
		return name;
	}

	/**
	 * The minimum score of this matrix.
	 * 
	 * @return minimum of the matrix.
	 */
	public short getMin() {
		return min;
	}

	/**
	 * The maximum score in this matrix.
	 * 
	 * @return maximum of the matrix.
	 */
	public short getMax() {
		return max;
	}

	/**
	 * Sets the description to the given value.
	 * 
	 * @param desc
	 *            a description. This doesn't have to start with '#'.
	 */
	public void setDescription(String desc) {
		this.description = desc;
	}

	/**
	 * Gives the alphabet used by this matrix.
	 * 
	 * @return the alphabet of this matrix.
	 */
	public FiniteAlphabet getAlphabet() {
		return alphabet;
	}

	/**
	 * Creates a <code>String</code> representation of this matrix.
	 * 
	 * @return a string representation of this matrix without the description.
	 */
	public String stringnifyMatrix() {
		int i = 0;
		StringBuffer matrixString = new StringBuffer();
		Symbol[] colSyms = new Symbol[this.colSymbols.keySet().size()];

		try {
			SymbolTokenization symtok = alphabet.getTokenization("default");
			matrixString.append("  ");
			Iterator<Symbol> colKeys = colSymbols.keySet().iterator();
			while (colKeys.hasNext()) {
				colSyms[i] = colKeys.next();
				matrixString.append(symtok.tokenizeSymbol(colSyms[i++])
						.toUpperCase());
				matrixString.append(' ');
			}
			matrixString.append(newLine);

			Iterator<Symbol> rowKeys = rowSymbols.keySet().iterator();
			while (rowKeys.hasNext()) {
				Symbol rowSym = rowKeys.next();
				matrixString
						.append(symtok.tokenizeSymbol(rowSym).toUpperCase());
				matrixString.append(' ');
				for (i = 0; i < colSyms.length; i++) {
					matrixString.append(getValueAt(rowSym, colSyms[i]));
					matrixString.append(' ');
				}
				matrixString.append(newLine);
			}
		} catch (BioException exc) {
			exc.printStackTrace();
		}
		return matrixString.toString();
	}

	/**
	 * Converts the description of the matrix to a String.
	 * 
	 * @return Gives a description with approximately 60 letters on every line
	 *         separated by <code>System.getProperty("line.separator")</code>.
	 *         Every line starts with <code>#</code>.
	 */
	public String stringnifyDescription() {
		StringBuffer desc = new StringBuffer(), line = new StringBuffer();
		line.append("# ");
		StringTokenizer st = new StringTokenizer(description, " ");
		while (st.hasMoreElements()) {
			line.append(st.nextElement().toString());
			line.append(' ');
			if (line.length() >= 60) {
				desc.append(line);
				desc.append(newLine);
				if (st.hasMoreElements()) {
					line = new StringBuffer();
					line.append("# ");
				}
			} else if (!st.hasMoreElements()) {
				desc.append(line);
				desc.append(newLine);
			}
		}
		return desc.toString();
	}

	/**
	 * Overrides the inherited method.
	 * 
	 * @return Gives a string representation of the SubstitutionMatrix. This is
	 *         a valid input for the constructor which needs a matrix string.
	 *         This String also contains the description of the matrix if there
	 *         is one.
	 */
	@Override
	public String toString() {
		StringBuffer desc = new StringBuffer(), line = new StringBuffer();
		line.append("# ");
		StringTokenizer st = new StringTokenizer(description);
		while (st.hasMoreElements()) {
			line.append(st.nextElement().toString());
			line.append(' ');
			if (line.length() >= 60) {
				desc.append(line);
				desc.append(newLine);
				if (st.hasMoreElements()) {
					line = new StringBuffer();
					line.append("# ");
				}
			} else if (!st.hasMoreElements()) {
				desc.append(line);
				desc.append(newLine);
			}
		}
		desc.append(stringnifyMatrix());
		return desc.toString();
	}

	/**
	 * Just to perform some test. It prints the matrix on the screen.
	 */
	public void printMatrix() {
		// Test output:
		Iterator<Symbol> rowKeys = rowSymbols.keySet().iterator();
		while (rowKeys.hasNext()) {
			Iterator<Symbol> colKeys = colSymbols.keySet().iterator();
			Symbol rowSym = rowKeys.next();
			System.out.print(rowSym.getName() + "\t");
			while (colKeys.hasNext()) {
				Symbol colSym = colKeys.next();
				int x = rowSymbols.get(rowSym).intValue();
				int y = colSymbols.get(colSym).intValue();
				System.out.print(colSym.getName() + " " + " " + x + " " + y
						+ " " + matrix[x][y] + "\t");
			}
			System.out.println(newLine);
		}
		System.out.println(toString());
	}

	/**
	 * With this method you can get a &ldquo;normalized&rdquo;
	 * <code>SubstitutionMatrix</code> object; however, since this
	 * implementation uses an short matrix, the normalized matrix will be scaled
	 * by ten. If you need values between zero and one, you have to divide every
	 * value returned by <code>getValueAt</code> by ten.
	 * 
	 * @return a new and normalized <code>SubstitutionMatrix</code> object given
	 *         by this substitution matrix. Because this uses an
	 *         <code>short</code> matrix, all values are scaled by 10.
	 * @throws BioException
	 * @throws IOException
	 * @throws NumberFormatException
	 */
	public SubstitutionMatrix normalizeMatrix() throws BioException,
			NumberFormatException, IOException {
		int i, j;
		short min = getMin(), newMax = Short.MIN_VALUE;
		short[][] mat = new short[matrix.length][matrix[matrix.length - 1].length];
		String name = getName() + "_normalized";
		String matString = stringnifyDescription() + "  ";
		FiniteAlphabet alphabet = getAlphabet();
		Map<Symbol, Integer> rowMap = this.rowSymbols;
		Map<Symbol, Integer> colMap = this.colSymbols;
		SymbolTokenization symtok = alphabet.getTokenization("default");

		for (i = 0; i < matrix.length; i++)
			for (j = 0; j < matrix[matrix.length - 1].length; j++) {
				mat[i][j] = (short) (matrix[i][j] - min);
				if (mat[i][j] > newMax)
					newMax = mat[i][j];
			}

		for (i = 0; i < mat.length; i++)
			for (j = 0; j < mat[mat.length - 1].length; j++)
				mat[i][j] = (short) (mat[i][j] * 10 / newMax);

		Object[] rows = rowSymbols.keySet().toArray();
		Object[] cols = colSymbols.keySet().toArray();
		for (i = 0; i < cols.length; i++)
			matString += symtok.tokenizeSymbol((Symbol) cols[i]) + " ";
		for (i = 0; i < rows.length; i++) {
			matString += newLine + symtok.tokenizeSymbol((Symbol) rows[i])
					+ " ";
			for (j = 0; j < cols.length; j++) {
				matString += mat[rowMap.get((Symbol) rows[i]).intValue()][colMap
						.get((Symbol) cols[j]).intValue()]
						+ " ";
			}
		}
		matString += newLine;
		return new SubstitutionMatrix(alphabet, matString, name);
	}

}
