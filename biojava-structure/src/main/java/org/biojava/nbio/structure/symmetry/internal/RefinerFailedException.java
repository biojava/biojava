package org.biojava.nbio.structure.symmetry.internal;

/**
 * Refinement of the symmetry alignment failed.
 * 
 * @author Aleix Lafita
 */
public class RefinerFailedException extends Exception {
	
	private static final long serialVersionUID = -3592155787060329421L;

	public RefinerFailedException() {
		super();
	}

	public RefinerFailedException(String message, Throwable cause) {
		super(message, cause);
	}

	public RefinerFailedException(String message) {
		super(message);
	}

	public RefinerFailedException(Throwable cause) {
		super(cause);
	}

}
