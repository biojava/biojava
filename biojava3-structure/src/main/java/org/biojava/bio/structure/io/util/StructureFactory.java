package org.biojava.bio.structure.io.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.cath.CathDatabase;
import org.biojava.bio.structure.cath.CathFactory;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;

/**
 * <p>
 * This class is intended as a global mechanism for retrieving {@link Structure}
 * objects given arbitrary identifiers. The set of recognized identifiers
 * <em>thus far</em> includes:
 * <ul>
 * <li>PDB Ids (ex: {@code 1hiv})</li>
 * <li>PDB Ids restricted to particular chains (ex: {@code 1hiv.A},
 * {@code 1hiv.a}, and {@code 1hiv.a,B})</li>
 * <li>PDB Ids restricted to particular explicit sequence ranges (ex:
 * {@code 1hiv.A_10-15}, {@code 1qdm.A_3S-99S}, and
 * {@code 1cph.B_4-6,B_11-30,A_1-9})</li>
 * <li>SCOP identifiers (ex: {@code d1hiva_})</li>
 * <li>CATH identifiers (ex: {@code 1cukA01})</li>
 * </ul>
 * More identifier types may be added in future versions, provided that they are
 * unambiguous. These identifiers may or may not be case-sensitive, though
 * case-insensitivity is preferred in cases where it is unambiguous.
 * </p>
 * <p>
 * This class implements a singleton pattern through {@link #getDefault()};
 * calling this method is preferred. This default uses "default" versions of
 * databases (such as {@link ScopFactory#DEFAULT_VERSION} and
 * {@link CathFactory#DEFAULT_VERSION}, and it is immutable by design. This
 * class also implements a multiton pattern through named key–value pairs; see
 * {@link #getInstance(String)} and
 * {@link #putInstance(String, StructureFactory)}. These latter functions are
 * intended only to facilitate unusual cases where distinct copies or versions
 * of databases are needed for different resources.
 * </p>
 * <p>
 * Note that this class is designed to live "above" its underlying databases.
 * This means that calling {@link ScopFactory#setScopDatabase(ScopDatabase)} and
 * {@link CathFactory#setCath(CathDatabase)}, for example, will have no effect.
 * </p>
 * 
 * @author dmyersturnbull
 * @since 3.0.8
 */
public class StructureFactory {

	private static StructureFactory defaultFactory;
	private static Map<String, StructureFactory> instances = new HashMap<String, StructureFactory>();

	/**
	 * Removes the StructureFactory with the specified identifier from the
	 * multiton. This method can be called in cases where an identifier must be
	 * removed to free memory for garbage collection.
	 */
	public static void destroyInstance(String factoryIdentifier) {
		instances.remove(factoryIdentifier);
	}

	/**
	 * Returns the default StructureFactory.
	 * Uses lazy initialization, so this <em>may take a while on the first call</em>.
	 * @return
	 */
	public static StructureFactory getDefault() {
		if (defaultFactory == null) {
			/*
			 * Use DEFAULT_VERSION because StructureFactory should be
			 * higher-level than either ScopFactory or CathFactory, so it should
			 * not be affected by setScop() or setCath()
			 */
			defaultFactory = new StructureFactory(new AtomCache(), ScopFactory.getSCOP(ScopFactory.DEFAULT_VERSION),
					CathFactory.getCathDatabase(CathFactory.DEFAULT_VERSION));
		}
		return defaultFactory;
	}

	/**
	 * From the multiton, returns the StructureFactory with name {@code factoryIdentifier}, or null if it doesn't exist.
	 * @see #putInstance(String, StructureFactory)
	 */
	public static StructureFactory getInstance(String factoryIdentifier) {
		return instances.get(factoryIdentifier);
	}

	/**
	 * Adds a named StructureFactory to the multiton. Where possible,
	 * <strong>calling {@link #getDefault()} is preferred</strong>.
	 * <em>If a StructureFactory with the specified identifier already exists, this method will do nothing.</em>
	 * 
	 * @param factoryIdentifier
	 *            A named identifier for the factory being added; this should be
	 *            sufficiently long as to avoid conflict
	 * @see #getInstance(String)
	 * @see #destroyInstance(String)
	 */
	public static void putInstance(String factoryIdentifier, StructureFactory newFactory) {
		if (!instances.containsKey(factoryIdentifier)) {
			instances.put(factoryIdentifier, newFactory);
		}
	}

	private final AtomCache cache;

	private final CathDatabase cath;

	private final ScopDatabase scop;

	/**
	 * Creates a new StructureFactory using the specified databases. Where
	 * possible, <strong>calling {@link #getDefault()} is preferred</strong>. A
	 * StructureFactory can be added statically to this class via
	 * {@link #putInstance(String, StructureFactory)}.
	 */
	public StructureFactory(AtomCache cache, ScopDatabase scop, CathDatabase cath) {
		this.cache = cache;
		this.scop = scop;
		this.cath = cath;
	}

	/**
	 * Returns underlying {@link AtomCache}.
	 */
	public AtomCache getAtomCache() {
		return cache;
	}

	/**
	 * Returns underlying {@link CathDatabase}.
	 */
	public CathDatabase getCath() {
		return cath;
	}

	/**
	 * Returns underlying {@link ScopDatabase}.
	 */
	public ScopDatabase getScop() {
		return scop;
	}

	/**
	 * Returns a {@link Structure} corresponding to the specified identifier.
	 * {@code identifier} can be of one of the following forms:
	 * <ul>
	 * <li>PDB Ids (ex: {@code 1hiv})</li>
	 * <li>PDB Ids restricted to particular chains (ex: {@code 1hiv.A},
	 * {@code 1hiv.a}, and {@code 1hiv.a,B})</li>
	 * <li>PDB Ids restricted to particular explicit sequence ranges (ex:
	 * {@code 1hiv.A_10-15}, {@code 1qdm.A_3S-99S}, and
	 * {@code 1cph.B_4-6,B_11-30,A_1-9})</li>
	 * <li>SCOP identifiers (ex: {@code d1hiva_})</li>
	 * <li>CATH identifiers (ex: {@code 1cukA01})</li>
	 * </ul>
	 */
	public Structure getStructure(String identifier) throws IOException, StructureException {
		StructureName theName = new StructureName(identifier);
		if (theName.isScopName()) {
			return cache.getStructureForDomain(identifier, scop);
		} else if (theName.isCathID()) {
			return cache.getStructureForCathDomain(theName, cath);
		} else {
			return cache.getStructure(identifier);
		}
	}

}
