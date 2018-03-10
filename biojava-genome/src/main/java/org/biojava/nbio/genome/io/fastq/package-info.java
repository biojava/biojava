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

/**
 * FASTQ and variants sequence format I/O.
 *
 * <p>
 * To read from an Illumina variant FASTQ sequence file:
 * <pre>
 * FastqReader reader = new IlluminaFastqReader();
 * for (Fastq fastq : reader.read(new File("illumina.fastq"))
 * {
 *   // ...
 * }
 * </pre>
 *
 * To write to a Sanger variant FASTQ sequence file:
 * <pre>
 * Collection&lt;Fastq&gt; fastq = ...;
 * SangerFastqWriter writer = new SangerFastqWriter();
 * writer.write(new File("sanger.fastq"), fastq);
 * </pre>
 *
 * For further documentation on the FASTQ sequence format,
 * its variants, and how they are handled in O|B|F projects,
 * see:
 *
 * <p>
 * 
 * <a href="http://dx.doi.org/10.1093/nar/gkp1137">The Sanger FASTQ file format for sequences
 * with quality scores, and the Solexa/Illumina FASTQ variants</a>
 * <p>
 * Peter J. A. Cock (Biopython), Christopher J. Fields (BioPerl), Naohisa Goto (BioRuby),
 * Michael L. Heuer (BioJava) and Peter M. Rice (EMBOSS).<br>
 * Nucleic Acids Research, <a href="http://dx.doi.org/10.1093/nar/gkp1137">doi:10.1093/nar/gkp1137</a>
 *
 * <p>
 * Moved from org.biojava.nbio.sequencing (biojava-sequencing module) in 5.0.0
 * 
 * @since 3.0.3
 * 
 */
package org.biojava.nbio.genome.io.fastq;
