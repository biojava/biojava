package org.biojava.nbio.structure.align.model;

/**
 * An exception during the handling or calculation of a strcuture alignment.
 * Created for the new Object Oriented {@link MultipleAlignment} Data Structure to store alignment results.
 * 
 * @author Aleix Lafita
 * 
 */
public class StructureAlignmentException extends Exception {

	private static final long serialVersionUID = -2886929241534584550L;

	public StructureAlignmentException() {
		super();
	}

	public StructureAlignmentException(String message, Throwable cause) {
		super(message, cause);
	}

	public StructureAlignmentException(String message) {
		super(message);
	}

	public StructureAlignmentException(Throwable cause) {
		super(cause);
	}

}