/**
 * 
 */
package org.biojava.bio.structure.align.benchmark;

/**
 * Interface for classes which parse alignments
 * 
 * Basically, the only requirement for a Parser is that it provide an iterator
 * to return one or more {@link MultipleAlignment}s
 * 
 * @author Spencer Bliven
 *
 */
public interface MultipleAlignmentParser extends Iterable<MultipleAlignment> {

}
