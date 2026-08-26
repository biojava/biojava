/*
 * BioJava development code
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
 */
package org.biojava.nbio.structure.ecod;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;

import org.biojava.nbio.structure.ecod.EcodInstallation.EcodParser;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;

/**
 * Checks that {@link EcodParser} reads every layout ECOD has distributed.
 * <p>
 * The columns have been renamed, reordered and added to several times, and the version
 * comment itself changed form at v294.1. Because the full distribution is 657 MB, none of
 * that was covered by a test that could run in reasonable time, and a format change went
 * unnoticed for months. These fixtures are taken verbatim from the real files, so the
 * contract is pinned in milliseconds rather than by a download.
 *
 * @author Amr ALHOSSARY
 * @since 7.3.0
 */
class EcodParserTest {

	/** develop204, list format 1.5: 15 columns, commented header, quoted names. */
	private static final String DEVELOP204 = String.join("\n",
			"#/data/ecod/database_versions/v204/ecod.develop204.domains.txt",
			"#ECOD version develop204",
			"#Domain list version 1.5",
			"#Grishin lab (http://prodata.swmed.edu/ecod)",
			"#uid\tecod_domain_id\tmanual_rep\tf_id\tpdb\tchain\tpdb_range\tseqid_range"
					+ "\tarch_name\tx_name\th_name\tt_name\tf_name\tasm_status\tligand",
			"002137905\te6b4nA1\tAUTO_NONREP\t1.1.1\t6b4n\tA\tA:1-99\tA:1-99\tbeta barrels"
					+ "\t\"cradle loop barrel\"\t\"RIFT-related\"\t\"acid protease\""
					+ "\tF_UNCLASSIFIED\tNOT_DOMAIN_ASSEMBLY\tCL,G53,NA");

	/** develop291, list format 1.6: 16 columns, f_id renamed t_id, unp_acc inserted at 9. */
	private static final String DEVELOP291 = String.join("\n",
			"#/data/ecod/database_versions/v291/ecod.develop291.domains.txt",
			"#ECOD version develop291",
			"#Domain list version 1.6",
			"#Grishin lab (http://prodata.swmed.edu/ecod)",
			"#uid\tecod_domain_id\tmanual_rep\tt_id\tpdb\tchain\tpdb_range\tseqid_range\tunp_acc"
					+ "\tarch_name\tx_name\th_name\tt_name\tf_name\tasm_status\tligand",
			"000000267\te1udzA1\tMANUAL_REP\t1.1.1\t1udz\tA\tA:203-381\tA:4-182\tP12345"
					+ "\tbeta barrels\t\"cradle loop barrel\"\t\"RIFT-related\"\t\"acid protease\""
					+ "\tF_UNCLASSIFIED\tNOT_DOMAIN_ASSEMBLY\tNO_LIGANDS_4A");

	private static final String V295_COLUMNS =
			"uid\tecod_domain_id\tmanual_rep\tf_id\tpdb\tchain\tpdb_range\tseqid_range"
			+ "\tarchitecture_name\tx_name\th_name\tt_name\tf_name\tassembly_id\tdomain_id_short"
			+ "\trange_count\tarch_manual\tx_manual\th_manual\tt_manual\tf_manual"
			+ "\tvalid_structure\tligand_binding\tligand_comp_ids\tligand_pdbnum";

	/** v295: 25 columns, uncommented header, True/False, empty assembly_id, moved ligands. */
	private static final String V295 = String.join("\n",
			"# ECOD Domain List",
			"# Version: v295",
			"# Generated: 2026-06-24 22:47:42",
			"# Ligand cutoff: 4.0 A (NO_LIGANDS_4A = no contact within cutoff)",
			"#",
			V295_COLUMNS,
			"0\te2nmzA1\tTrue\t1.1.1.3\t2nmz\tA\tA:1-99\tA:1-99\tbeta barrels\tcradle loop barrel"
					+ "\tRIFT-related\tacid protease\tRVP\t\t\t1\tFalse\tFalse\tFalse\tFalse\tTrue"
					+ "\tTrue\tTrue\tROC,SO4\tA:601,A:602,B:401",
			// the last column is empty on four rows in five, so split() must keep it
			"3\te2rspA1\tTrue\t1.1.1.3\t2rsp\tA\tA:1-124\tA:1-124\tbeta barrels\tcradle loop barrel"
					+ "\tRIFT-related\tacid protease\tRVP\t\t\t1\tFalse\tFalse\tFalse\tFalse\tTrue"
					+ "\tTrue\tFalse\tNO_LIGANDS_4A\t",
			// a domain classified from an AlphaFold model: no PDB entry, so no EcodDomain
			"3163557\tP44140_F1_nD2\tFalse\t2004.1.1.123\t\t\t131-315\t131-315\talpha bundles"
					+ "\tsomething\tsomething else\ta third thing\t\t\t\t1\tFalse\tFalse\tFalse"
					+ "\tFalse\tTrue\tTrue\tFalse\tNO_LIGANDS_4A\t");

	private static List<EcodDomain> parse(String contents) throws IOException {
		return new EcodParser(new StringReader(contents)).getDomains();
	}

	private static String version(String contents) throws IOException {
		return new EcodParser(new StringReader(contents)).getVersion();
	}

	@Nested
	class Version {
		@Test
		void oldHeaderForm() throws IOException {
			assertEquals("develop204", version(DEVELOP204));
			assertEquals("develop291", version(DEVELOP291));
		}

		@Test
		void newHeaderForm() throws IOException {
			assertEquals("v295", version(V295));
			assertEquals("v294.1", version("# ECOD Domain List\n# Version: v294.1\n"));
		}

		@Test
		void listFormatVersionIsNotTheEcodVersion() throws IOException {
			// "#Domain list version 1.5" describes the columns, not the release
			assertNull(version("#Grishin lab\n#Domain list version 1.5\n"));
		}

		@Test
		void absentVersionIsNull() throws IOException {
			assertNull(version("#Grishin lab (http://prodata.swmed.edu/ecod)\n"));
		}
	}

	@Nested
	class OldFormats {
		@Test
		void listFormat15() throws IOException {
			List<EcodDomain> domains = parse(DEVELOP204);
			assertEquals(1, domains.size());
			EcodDomain d = domains.get(0);
			assertEquals(Long.valueOf(2137905), d.getUid());
			assertEquals("e6b4nA1", d.getDomainId());
			assertEquals(Boolean.FALSE, d.getManual());
			assertEquals(Integer.valueOf(1), d.getXGroup());
			assertEquals(Integer.valueOf(1), d.getHGroup());
			assertEquals(Integer.valueOf(1), d.getTGroup());
			assertNull(d.getFGroup());
			assertEquals("6B4N", d.getPdbId().getId());
			assertEquals("A", d.getChainId());
			assertEquals("A:1-99", d.getRange());
			assertEquals("A:1-99", d.getSeqIdRange());
			assertEquals("beta barrels", d.getArchitectureName());
			// quotes were stripped up to develop292
			assertEquals("cradle loop barrel", d.getXGroupName());
			assertEquals("RIFT-related", d.getHGroupName());
			assertEquals("acid protease", d.getTGroupName());
			assertEquals("F_UNCLASSIFIED", d.getFGroupName());
			assertEquals(Long.valueOf(2137905), d.getAssemblyId());
			assertEquals(new LinkedHashSet<>(Arrays.asList("CL", "G53", "NA")), d.getLigands());
		}

		/**
		 * develop291 inserted unp_acc before arch_name. Read positionally, every field from
		 * there on shifts by one and the domain is silently mangled or dropped.
		 */
		@Test
		void listFormat16InsertsUnpAcc() throws IOException {
			List<EcodDomain> domains = parse(DEVELOP291);
			assertEquals(1, domains.size());
			EcodDomain d = domains.get(0);
			assertEquals(Boolean.TRUE, d.getManual());
			assertEquals("1UDZ", d.getPdbId().getId());
			assertEquals("A:4-182", d.getSeqIdRange());
			assertEquals("beta barrels", d.getArchitectureName());
			assertEquals("acid protease", d.getTGroupName());
			assertEquals("F_UNCLASSIFIED", d.getFGroupName());
			assertEquals(Long.valueOf(267), d.getAssemblyId());
			assertEquals(Collections.emptySet(), d.getLigands());
		}

		/**
		 * Headers were only added in develop101. Older files are still read by position.
		 */
		@Test
		void thirteenColumnsWithoutAHeader() throws IOException {
			List<EcodDomain> domains = parse(String.join("\n",
					"#ECOD version develop45",
					"000000001\te1udzA1\t1.1.1\t1udz\tA\tA:203-381\tbeta barrels"
							+ "\t\"cradle loop barrel\"\t\"RIFT-related\"\t\"acid protease\""
							+ "\tF_UNCLASSIFIED\tNOT_DOMAIN_ASSEMBLY\tNO_LIGANDS_4A"));
			assertEquals(1, domains.size());
			EcodDomain d = domains.get(0);
			assertNull(d.getManual(), "no manual_rep column before list format 1.1");
			assertNull(d.getSeqIdRange(), "no seqid_range column before list format 1.4");
			assertEquals("1UDZ", d.getPdbId().getId());
			assertEquals("acid protease", d.getTGroupName());
		}
	}

	@Nested
	class NewFormat {
		@Test
		void twentyFiveColumns() throws IOException {
			List<EcodDomain> domains = parse(V295);
			// two PDB domains; the AlphaFold-derived row cannot be an EcodDomain
			assertEquals(2, domains.size());

			EcodDomain d = domains.get(0);
			assertEquals(Long.valueOf(0), d.getUid());
			assertEquals("e2nmzA1", d.getDomainId());
			assertEquals(Boolean.TRUE, d.getManual(), "manual_rep is now True/False");
			assertEquals(Integer.valueOf(1), d.getXGroup());
			assertEquals(Integer.valueOf(1), d.getHGroup());
			assertEquals(Integer.valueOf(1), d.getTGroup());
			assertEquals(Integer.valueOf(3), d.getFGroup(), "f_id now carries a fourth level");
			assertEquals("2NMZ", d.getPdbId().getId());
			assertEquals("A", d.getChainId());
			assertEquals("A:1-99", d.getRange());
			assertEquals("A:1-99", d.getSeqIdRange());
			assertEquals("beta barrels", d.getArchitectureName());
			assertEquals("cradle loop barrel", d.getXGroupName());
			assertEquals("RIFT-related", d.getHGroupName());
			assertEquals("acid protease", d.getTGroupName());
			assertEquals("RVP", d.getFGroupName());
			// assembly_id is empty on every row of v294.1 and later, which means the same
			// as the NOT_DOMAIN_ASSEMBLY of earlier versions
			assertEquals(Long.valueOf(0), d.getAssemblyId());
			assertEquals(new LinkedHashSet<>(Arrays.asList("ROC", "SO4")), d.getLigands(),
					"the ligand list moved to ligand_comp_ids");
		}

		/**
		 * ligand_pdbnum, the last column, is empty on four rows in five. String.split
		 * discards trailing empty fields unless asked not to, which would make those rows
		 * look one column short.
		 */
		@Test
		void rowEndingInAnEmptyColumn() throws IOException {
			EcodDomain d = parse(V295).get(1);
			assertEquals("e2rspA1", d.getDomainId());
			assertEquals("2RSP", d.getPdbId().getId());
			assertEquals(Collections.emptySet(), d.getLigands());
		}

		@Test
		void twentyThreeColumnsOfV2941() throws IOException {
			List<EcodDomain> domains = parse(String.join("\n",
					"# ECOD Domain List",
					"# Version: v294.1",
					"#",
					"uid\tecod_domain_id\tmanual_rep\tf_id\tpdb\tchain\tpdb_range\tseqid_range"
							+ "\tarchitecture_name\tx_name\th_name\tt_name\tf_name\tassembly_id"
							+ "\tdomain_id_short\trange_count\tarch_manual\tx_manual\th_manual"
							+ "\tt_manual\tf_manual\tvalid_structure\tligand_binding",
					"1\te1hvcA1\tFalse\t1.1.1.3\t1hvc\tA\tA:1B-99A\tA:1-203\tbeta barrels"
							+ "\tcradle loop barrel\tRIFT-related\tacid protease\tRVP\t\t\t\tFalse"
							+ "\tFalse\tFalse\tFalse\tFalse\tTrue\tFalse"));
			assertEquals(1, domains.size());
			EcodDomain d = domains.get(0);
			assertEquals("1HVC", d.getPdbId().getId());
			assertEquals("A:1B-99A", d.getRange());
			assertEquals("RVP", d.getFGroupName());
			// there is no ligand column at all in v294.1
			assertEquals(Collections.emptySet(), d.getLigands());
		}

		@Test
		void columnHeaderIsNotADomain() throws IOException {
			// v294.1 stopped commenting the column names out, so they arrive looking like data
			for (EcodDomain d : parse(V295)) {
				assertFalse("uid".equals(d.getDomainId()));
			}
		}

		/**
		 * An empty f_name is not the same as F_UNCLASSIFIED: f_id still classifies the
		 * domain to four levels, so the empty value is left as it is rather than translated.
		 */
		@Test
		void emptyFGroupNameIsLeftAlone() throws IOException {
			List<EcodDomain> domains = parse(String.join("\n",
					"# Version: v295",
					V295_COLUMNS,
					"7\te4fivA1\tTrue\t1.1.1.3\t4fiv\tA\tA:4-116\tA:1-113\tbeta barrels"
							+ "\tcradle loop barrel\tRIFT-related\tacid protease\t\t\t\t1\tFalse"
							+ "\tFalse\tFalse\tFalse\tTrue\tTrue\tTrue\tLP1\tA:201"));
			assertEquals(1, domains.size());
			assertEquals("", domains.get(0).getFGroupName());
		}
	}

	@Nested
	class Robustness {
		@Test
		void unparseableLinesAreSkippedNotFatal() throws IOException {
			List<EcodDomain> domains = parse(String.join("\n",
					"# Version: v295",
					V295_COLUMNS,
					"not-a-number\tefoo\tTrue\t1.1.1.3\tfoo1\tA\tA:1-9\tA:1-9\ta\tb\tc\td\te"
							+ "\t\t\t1\tFalse\tFalse\tFalse\tFalse\tTrue\tTrue\tFalse\t\t",
					"0\te2nmzA1\tTrue\t1.1.1.3\t2nmz\tA\tA:1-99\tA:1-99\tbeta barrels"
							+ "\tcradle loop barrel\tRIFT-related\tacid protease\tRVP\t\t\t1"
							+ "\tFalse\tFalse\tFalse\tFalse\tTrue\tTrue\tTrue\tROC,SO4\tA:601"));
			assertEquals(1, domains.size(), "the good line is still read");
			assertEquals("e2nmzA1", domains.get(0).getDomainId());
		}

		@Test
		void shortLineIsSkipped() throws IOException {
			assertTrue(parse(String.join("\n",
					"# Version: v295",
					V295_COLUMNS,
					"0\te2nmzA1\tTrue")).isEmpty());
		}

		@Test
		void emptyFileYieldsNoDomains() throws IOException {
			assertTrue(parse("").isEmpty());
		}
	}
}
