package org.biojava.nbio.structure.io.cif;

/**
 * Defines a rather generic interface which allows to populate some data structure with data parsed from a CIF file.
 * @param <S> the type of container an implementing class will return
 * @author Sebastian Bittrich
 * @since 5.3.0
 */
public interface CifFileConsumer<S> {
    /**
     * Setup routine which initializes a new container.
     */
    void prepare();

    /**
     * Ultimate setup which can include steps which require several categories to be available and integrate them into
     * the final container.
     */
    void finish();

    /**
     * Retrieve the created container representing a CIF file.
     * @return all desired information wrapped as object of type <code>S</code>
     */
    S getContainer();
}
