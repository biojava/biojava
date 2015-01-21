package org.biojava.bio.structure.scop;

/**
 * Indicates that an I/O error occurred with SCOP lazy initialization.
 */
public class ScopIOException extends RuntimeException {

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

	public ScopIOException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}
}
