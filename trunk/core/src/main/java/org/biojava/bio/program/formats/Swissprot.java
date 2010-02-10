package org.biojava.bio.program.formats;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.RegexFieldFinder;
import org.biojava.bio.program.tagvalue.RegexSplitter;
import org.biojava.bio.program.tagvalue.SimpleTagValueWrapper;
import org.biojava.bio.program.tagvalue.TagDelegator;
import org.biojava.bio.program.tagvalue.TagValueContext;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ParserException;
import org.biojava.utils.lsid.LifeScienceIdentifier;

public class Swissprot
implements Format {
  private static final AnnotationType ANNO_TYPE;
  //private static final LineSplitParser PARSER;
  private static final LifeScienceIdentifier LSID;

  static {
    LSID = LifeScienceIdentifier.valueOf("open-bio.org", "format", "swissprot");

    Location NONE = CardinalityConstraint.NONE;
    Location ANY = CardinalityConstraint.ANY;
    Location ONE = CardinalityConstraint.ONE;
    Location ONE_OR_MORE = CardinalityConstraint.ONE_OR_MORE;

    //PARSER = new LineSplitParser(LineSplitParser.EMBL);

    PropertyConstraint c_string = new PropertyConstraint.ByClass(String.class);

    AnnotationType FT = new AnnotationType.Impl();
    FT.setDefaultConstraints(PropertyConstraint.ANY, ANY); // fix this
    PropertyConstraint c_ft = new PropertyConstraint.ByAnnotationType(FT);

    ANNO_TYPE = new AnnotationType.Impl();
    ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
    ANNO_TYPE.setConstraints("ID", c_string, ONE);
    ANNO_TYPE.setConstraints("TYPE", c_string, ONE);
    ANNO_TYPE.setConstraints("MOLECULE", c_string, ONE);
    ANNO_TYPE.setConstraints("LENGTH", c_string, ONE);
    ANNO_TYPE.setConstraints("AC", c_string, ONE_OR_MORE);
    ANNO_TYPE.setConstraints("DT", c_string, ANY);
    ANNO_TYPE.setConstraints("KW", c_string, ANY);
    ANNO_TYPE.setConstraints("OS", c_string, ONE);
    ANNO_TYPE.setConstraints("OC", c_string, ANY);
    ANNO_TYPE.setConstraints("DE", c_string, ANY);
    ANNO_TYPE.setConstraints("GN", c_string, ANY);
    ANNO_TYPE.setConstraints("OS", c_string, ANY);
    ANNO_TYPE.setConstraints("OG", c_string, ANY);
    ANNO_TYPE.setConstraints("OC", c_string, ANY);
    ANNO_TYPE.setConstraints("OX", c_string, ANY);
    ANNO_TYPE.setConstraints("RN", c_string, ANY);
    ANNO_TYPE.setConstraints("RP", c_string, ANY);
    ANNO_TYPE.setConstraints("RC", c_string, ANY);
    ANNO_TYPE.setConstraints("RX", c_string, ANY);
    ANNO_TYPE.setConstraints("RA", c_string, ANY);
    ANNO_TYPE.setConstraints("RT", c_string, ANY);
    ANNO_TYPE.setConstraints("RL", c_string, ANY);
    ANNO_TYPE.setConstraints("CC", c_string, ANY);
    ANNO_TYPE.setConstraints("DR", c_string, ANY);
    ANNO_TYPE.setConstraints("KW", c_string, ANY);
    ANNO_TYPE.setConstraints("FT", c_ft, ANY);
    ANNO_TYPE.setConstraints("SQ", c_string, ANY);
    ANNO_TYPE.setConstraints("", c_string, ANY);
  }

  public ParserListener getParserListener(TagValueListener listener) {
    RegexSplitter semiColonSplitter = new RegexSplitter(
      Pattern.compile("(\\w+)[;.]"),
      1
    );
    ValueChanger semiColonChanger = new ValueChanger(listener);
    semiColonChanger.setDefaultSplitter(semiColonSplitter);

    LineSplitParser ftParser = new LineSplitParser();
    ftParser.setSplitOffset(29);
    ftParser.setTrimTag(true);
    ftParser.setTrimValue(true);
    ftParser.setContinueOnEmptyTag(true);
    ftParser.setMergeSameTag(false);

    TagValueListener ftListener = new SPFeatureTableListener(listener);

    LineSplitParser lsp = LineSplitParser.EMBL;
    TagDelegator td = new TagDelegator(listener);

    td.setListener("ID", new RegexFieldFinder(
      listener,
      Pattern.compile("(\\w+)\\s+(\\w+);\\s+(\\w+);\\s+(\\d+)"),
      new String[] { "ID", "TYPE", "MOLECULE", "LENGTH" },
      true
    ));
    td.setListener("AC", semiColonChanger);
    td.setListener("KW", semiColonChanger);
    td.setListener("OC", semiColonChanger);
    td.setListener("RC", semiColonChanger);
    td.setListener("RX", semiColonChanger);
    td.setParserListener("FT", ftParser, ftListener);

    return new ParserListener(lsp, td);
  }


  public AnnotationType getType() {
    return ANNO_TYPE;
  }

  public LifeScienceIdentifier getLSID() {
    return LSID;
  }

  private static class SPFeatureTableListener
  extends SimpleTagValueWrapper {
    private Pattern pat = Pattern.compile("(\\w+)\\s+((<?\\d+)|(?))\\s+((>?\\d+)|(\\?))");
    private int depth = 0;
    private Object tag;

    public SPFeatureTableListener(TagValueListener delegate) {
      super(delegate);
    }

    public void startRecord()
    throws ParserException {
      depth++;
      super.startRecord();
    }

    public void endRecord()
    throws ParserException {
      super.endRecord();
      depth--;
    }

    public void startTag(Object tag)
    throws ParserException {
      if(depth == 1) {
        this.tag = tag;
      } else {
        super.startTag(tag);
      }
    }

    public void endTag(Object tag)
    throws ParserException {
      if(depth == 1) {
        // do we need something here?
      }

      super.endTag();
    }

    public void value(TagValueContext ctxt, Object val)
    throws ParserException {
      if(depth == 1) {
        if(tag != null) {
          try {
            Matcher m = pat.matcher(tag.toString());
            m.find();

            super.startTag("TYPE");
            super.value(ctxt, m.group(1));
            super.endTag();

            super.startTag("START");
            super.value(ctxt, m.group(2));
            super.endTag();

            super.startTag("END");
            super.value(ctxt, m.group(3));
            super.endTag();

            super.startTag("DESCRIPTION");
            super.value(ctxt, val);

            tag = null;
          } catch (IllegalStateException ise) {
            throw new ParserException("Couldn't match: " + pat.pattern() + " " + tag, ise);
          }
        } else {
          super.value(ctxt, val);
        }
      } else {
        super.value(ctxt, val);
      }
    }
  }
}

