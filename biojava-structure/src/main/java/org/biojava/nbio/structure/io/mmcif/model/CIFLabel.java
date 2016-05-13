package org.biojava.nbio.structure.io.mmcif.model;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Annotation indicating that a specific field of a bean should be mapped to
 * a different label
 * @author Spencer Bliven
 *
 */
@Target(value=ElementType.FIELD)
@Retention(value=RetentionPolicy.RUNTIME)
public @interface CIFLabel {
	String label();
}
