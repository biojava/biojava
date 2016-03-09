/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.ontology;

import org.biojava.nbio.ontology.io.TabDelimParser;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Tools for manipulating ontologies.
 *
 * @author Matthew Pocock
 */
public final class OntoTools {

	private static final Logger logger = LoggerFactory.getLogger(OntoTools.class);

	private static final Ontology CORE_ONTOLOGY;
	private static final OntologyFactory DEFAULT_FACTORY;
	private static final IntegerOntology CORE_INTEGER;
	//private static final Ontology CORE_STRING;

	// public static final Term TYPE;
	public static final Term RELATION;
	public static final Term ANY;
	public static final Term NONE;
	public static final Term IS_A;
	public static final Term PART_OF;
	// public static final Term SUB_TYPE_OF;
	// public static final Term INSTANCE_OF;
	// public static final Term DOMAIN;
	// public static final Term CO_DOMAIN;
	// public static final Term HAS_DOMAIN;
	// public static final Term HAS_CO_DOMAIN;

	// public static final Term BOOLEAN;
	// public static final Term TRUE;
	// public static final Term FALSE;
	// public static final Term PREDICATE;

	// public static final Term AND;
	// public static final Term OR;
	// public static final Term XOR;
	// public static final Term EQUAL;
	// public static final Term NOT_EQUAL;
	// public static final Term IMPLIES;

	public static final Term REFLEXIVE;
	public static final Term SYMMETRIC;
	public static final Term TRANSITIVE;
	public static final Term EQUIVALENCE;
	public static final Term PARTIAL_ORDER;

	static {
		DEFAULT_FACTORY = new OntologyFactory() {
			@Override
			public Ontology createOntology(String name, String desc)
			throws OntologyException {
				return new Ontology.Impl(name, desc);
			}
		};

		try {
			BufferedReader reader = new BufferedReader(
				new InputStreamReader(
					OntoTools.class.getResourceAsStream(
						"/ontology/core.onto"
					)
				)
			);

			CORE_INTEGER = new IntegerOntology();
			CORE_ONTOLOGY = new TabDelimParser().parse(
							reader,
							DEFAULT_FACTORY
			);

			// TYPE = CORE_ONTOLOGY.getTerm("type");
			RELATION = CORE_ONTOLOGY.getTerm("relation");
			ANY = CORE_ONTOLOGY.getTerm("any");
			NONE = CORE_ONTOLOGY.getTerm("none");
			IS_A = CORE_ONTOLOGY.getTerm("is-a");
			PART_OF = CORE_ONTOLOGY.getTerm("part-of");

			// SUB_TYPE_OF = CORE_ONTOLOGY.getTerm("sub_type_of");
			// INSTANCE_OF = CORE_ONTOLOGY.getTerm("instance_of");
			// DOMAIN = CORE_ONTOLOGY.getTerm("domain");
			// CO_DOMAIN = CORE_ONTOLOGY.getTerm("co_domain");
			// HAS_DOMAIN = CORE_ONTOLOGY.getTerm("has_domain");
			// HAS_CO_DOMAIN = CORE_ONTOLOGY.getTerm("has_co_domain");

			// BOOLEAN = CORE_ONTOLOGY.getTerm("boolean");
			// TRUE = CORE_ONTOLOGY.getTerm("true");
			// FALSE = CORE_ONTOLOGY.getTerm("false");
			// PREDICATE = CORE_ONTOLOGY.getTerm("predicate");

			// AND = CORE_ONTOLOGY.getTerm("and");
			// OR = CORE_ONTOLOGY.getTerm("or");
			// XOR = CORE_ONTOLOGY.getTerm("xor");
			// EQUAL = CORE_ONTOLOGY.getTerm("equal");
			// NOT_EQUAL = CORE_ONTOLOGY.getTerm("not_equal");
			// IMPLIES = CORE_ONTOLOGY.getTerm("implies");

			REFLEXIVE = CORE_ONTOLOGY.getTerm("reflexive");
			EQUIVALENCE = CORE_ONTOLOGY.getTerm("equivalence");
			SYMMETRIC = CORE_ONTOLOGY.getTerm("symmetric");
			TRANSITIVE = CORE_ONTOLOGY.getTerm("transitive");
			PARTIAL_ORDER = CORE_ONTOLOGY.getTerm("partial-order");

		} catch (Exception e) {
			logger.error("Exception: ", e);
			throw new RuntimeException("Could not initialize OntoTools", e);
		}
	}


	private OntoTools() {}

	/**
	 * Get the Ontology that defines our core "central dogma".
	 *
	 * <p>This contains definitions that we have to have, such as <code>any</code>,
	 * <code>predicate</code>, <code>is-a</code> and <code>transient</code>. These
	 * are our axioms, upon which the default interpreters build.</p>
	 *
	 * @return the "core" Ontology
	 */
	public static Ontology getCoreOntology() {
		return CORE_ONTOLOGY;
	}

	/**
	 * Get the Ontology that defines integers.
	 *
	 * <p>This contains a term for each and every integer. I haven't decided yet
	 * if it contains terms for arithmatic.</p>
	 *
	 * @return the integer Ontology
	 */
	public static IntegerOntology getIntegerOntology() {
		return CORE_INTEGER;
	}

	public static OntologyFactory getDefaultFactory() {
		return DEFAULT_FACTORY;
	}
}
