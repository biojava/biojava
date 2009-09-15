package org.biojava.bio.program.formats;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.lsid.LifeScienceIdentifier;

public class Embl
implements Format {
  private static final AnnotationType ANNO_TYPE;
  private static final LineSplitParser PARSER;
  private static final LifeScienceIdentifier LSID;

  static {
    LSID = LifeScienceIdentifier.valueOf("open-bio.org", "format", "embl");

    Location NONE = CardinalityConstraint.NONE;
    Location ANY = CardinalityConstraint.ANY;
    Location ONE = CardinalityConstraint.ONE;
    Location ONE_OR_MORE = CardinalityConstraint.ONE_OR_MORE;

    PARSER = new LineSplitParser(LineSplitParser.EMBL);

    PropertyConstraint c_string = new PropertyConstraint.ByClass(String.class);

    ANNO_TYPE = new AnnotationType.Impl();
    ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
    ANNO_TYPE.setConstraints("ID", c_string, ONE);
    ANNO_TYPE.setConstraints("AC", c_string, ONE_OR_MORE);
    ANNO_TYPE.setConstraints("SV", c_string, ONE);
    ANNO_TYPE.setConstraints("DT", c_string, ANY);
  }

  public LifeScienceIdentifier getLSID() {
    return LSID;
  }

  public AnnotationType getType() {
    return ANNO_TYPE;
  }

  public ParserListener getParserListener(TagValueListener listener) {
    return new ParserListener(PARSER, listener);
  }
}
