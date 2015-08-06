package org.biojava.nbio.structure.symmetry.internal;

/**
 * Refinement of the self-alignment failed.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 * 
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
