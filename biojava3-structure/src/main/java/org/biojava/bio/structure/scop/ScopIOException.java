package org.biojava.bio.structure.scop;

/**
 * Indicates that an I/O error occurred with SCOP lazy initialization.
 */
public class ScopIOException extends RuntimeException {

	private static final long serialVersionUID = 1L;

	public ScopIOException() {
	}

	public ScopIOException(String message) {
		super(message);
	}

	public ScopIOException(String message, Throwable cause) {
		super(message, cause);
	}

	public ScopIOException(Throwable cause) {
		super(cause);
	}

}
