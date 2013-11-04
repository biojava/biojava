package org.biojava.bio.structure.scop;

/**
 * Classes which implement ScopDatabase in a way which allows them to serve queries
 * without accessing the internet should implement this interface instead. An
 * initial file download is acceptable.
 * 
 * ScopFactory utilizes this distinction to optimize when many queries are expected
 * or remote calls would be otherwise undesireable.
 * @author Spencer Bliven
 *
 */
public interface LocalScopDatabase extends ScopDatabase {}
