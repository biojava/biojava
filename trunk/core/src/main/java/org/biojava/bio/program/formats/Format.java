package org.biojava.bio.program.formats;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.utils.lsid.LifeScienceIdentifier;

/**
 * A file format supported by the tag-value event-based parsing system.
 *
 * <p>Format should be implemented to provide pre-canned access to common
 * formats such as Embl, genbank, swissprot, enzyme etc. so that people do not
 * need to work out which events should be produced by a given file format.
 * It is expcected that implementations of Format will publish meta-data
 * about what tags are associated with which values.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public interface Format {
  /**
   * Retrieve a ParserListener pair for the format that will pass all events
   * on to a listener.
   *
   * Call this method to get a working parser that can be fed into a tag-value
   * parsing pipeline.
   *
   * This method may well be called many times during the lifetime of an
   * applications.. You should make this threadsafe. To avoid buring too much
   * memory, and to facilitate the comparrison of object by the == operator,
   * it is usefull to shair as much data as possible between the parsers and
   * handlers returned by this method.
   *
   * @param listener  a TagValueListener that all events should be passed onto
   * @return a ParserListener for the format
   */
  ParserListener getParserListener(TagValueListener listener);

  /**
   * Get the AnnotationType that constrains the events that will be fired.
   *
   * If you feed the events from the ParserListener into somethign that builds
   * Annotation bundles, this is the AnnotationType that those bundles will
   * conform to.
   *
   * In the cases where the events have been sensibly crafted, it will be
   * possible to introspect a great deal about the parsing events from this
   * AnntoationType. Use it to dynamicaly bind events to object models, generate
   * gui componets, and to work out which formats contain cross-refferenceable
   * information.
   *
   * It is polite to return a full and constrained description of the types of
   * oevents that may be generated, how many of them could come (cardinality)
   * and what types of values will be associated with them. The use of
   * OmtologyTerm instances as property names is encouraged.
   *
   * @return an AnnotationType representingchema for the events
   */
  AnnotationType getType();

  /**
   * Retrieve the LSID associated with this format.
   *
   * <p>The OBDA recomends taht file formats have identifiers assopciated with
   * them. This allows the format to be specified unambiguously across
   * different projects and groups. Idealy, a format LSID should conform to
   * the odda <a href="http://cvs.open-bio.org/cgi-bin/viewcvs/viewcvs.cgi/obda-specs/registry/lsid_for_dbformats.txt?rev=HEAD&cvsroot=obf-common&content-type=text/vnd.viewcvs-markup">formats specification</a>.</p>
   */
  LifeScienceIdentifier getLSID();
}

