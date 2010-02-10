package org.biojava.bio.program.tagvalue;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.CollectionConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.utils.ParserException;

/**
 * This is intended as a repository for tag-value and AnnotationType information
 * about common file formats. Each format should have an annotaiton type
 * defined as &lt;FormatName&gt;_TYPE and a method
 * create&lt;FormatName&gt;ParserListener(ParserListener listener) that together
 * give you everything needed to parse and represent the format.
 *
 * @author Matthew Pocock
 */
public class Formats {
  public static final AnnotationType EMBL_TYPE;
  public static final AnnotationType EMBL_GENBANK_FEATURE_TABLE_TYPE;
  public static final AnnotationType SWISSPROT_TYPE;

  static {
    PropertyConstraint prop_string = new PropertyConstraint.ByClass(String.class);
    CollectionConstraint prop_stringList = new CollectionConstraint.AllValuesIn(
      prop_string,
      CardinalityConstraint.ANY
    );

    // feature table strucure - shared by embl & genbank
    EMBL_GENBANK_FEATURE_TABLE_TYPE = new AnnotationType.Impl();
    EMBL_GENBANK_FEATURE_TABLE_TYPE.setDefaultConstraint(prop_stringList);
    PropertyConstraint prop_featureTable = new PropertyConstraint.ByAnnotationType(EMBL_GENBANK_FEATURE_TABLE_TYPE);

    // embl top-level
    EMBL_TYPE = new AnnotationType.Impl();
    EMBL_TYPE.setDefaultConstraint(prop_stringList);
    EMBL_TYPE.setConstraints("FT", prop_featureTable, CardinalityConstraint.ZERO_OR_ONE);

    // swissprot top-level
    SWISSPROT_TYPE = new AnnotationType.Impl();
    SWISSPROT_TYPE.setDefaultConstraint(prop_stringList);
  }

  public static final ParserListener createEmblParserListener(TagValueListener listener) {
    RegexSplitter semiColonSplitter = new RegexSplitter(
      Pattern.compile("(\\w+)[;.]"),
      1
    );
    ValueChanger semiColonChanger = new ValueChanger(listener);
    semiColonChanger.setDefaultSplitter(semiColonSplitter);


    LineSplitParser lsp = LineSplitParser.EMBL;

    TagDelegator td = new TagDelegator(listener);

    LineSplitParser ftParser = new LineSplitParser();
    ftParser.setSplitOffset(15);
    ftParser.setTrimTag(true);
    ftParser.setTrimValue(true);
    ftParser.setContinueOnEmptyTag(true);
    ftParser.setMergeSameTag(false);

    TagValueListener ftListener = new FeatureTableListener(listener);

    td.setParserListener("FT", ftParser, ftListener);
    td.setListener("ID", new RegexFieldFinder(
      listener,
      Pattern.compile("(\\w+)\\s+(\\w+);\\s+(.*?);\\s+(\\w+);\\s+(\\d+)\\s+BP\\."),
      new String[] { "ID", "TYPE", "MOLECULE", "DIVISION", "SIZE" },
      true
    ));
    td.setListener("AC", semiColonChanger);
    td.setListener("KW", semiColonChanger);
    td.setListener("OC", semiColonChanger);

    return new ParserListener(lsp, td);
  }

  public static final ParserListener createSwissprotParserListener(TagValueListener listener) {
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

  private static class FeatureTableListener
  extends SimpleTagValueWrapper {
    private TagValueParser featurePropertyParser = new FeaturePropertyParser();
    private int depth = 0;
    
    private boolean inLocation;

    public FeatureTableListener() {
      super();
    }

    public FeatureTableListener(TagValueListener delegate) {
      super(delegate);
    }

    public void startRecord()
    throws ParserException  {
      inLocation = false;

      super.startRecord();
    }

    public void endRecord()
    throws ParserException {
      if(inLocation) {
        super.endTag();
      }

      super.endRecord();
    }

    public void startTag(Object tag)
    throws ParserException {
      super.startTag(tag);

      if(depth == 0) {
        super.startRecord();
      }

      depth++;
    }

    public void endTag()
    throws ParserException {
      depth--;

      if(depth == 0) {
        super.endRecord();
      }

      super.endTag();
    }

    public void value(TagValueContext tvc, Object value)
    throws ParserException {
      String line = (String) value;
      if(line.startsWith("/")) {
        if(inLocation) {
          super.endTag();
          inLocation = false;
        }
        tvc.pushParser(featurePropertyParser, new TopRecordDropper(getDelegate()));
      } else {
        if(!inLocation) {
          super.startTag("LOCATION");
          inLocation = true;
        }
        super.value(tvc, value);
      }
    }
  }

  private static class FeaturePropertyParser
  implements TagValueParser {
    public TagValue parse(Object value)
    throws ParserException  {
      String line = (String) value;
      if(line.startsWith("/")) {
        int eq = line.indexOf("=");
        if(eq < 0) {
          return new TagValue(line.substring(1), "", true);
        } else {
          String ourTag = line.substring(1, eq);
          String ourValue = line.substring(eq + 1);
          return new TagValue(ourTag, ourValue, true);
        }
      } else {
        return new TagValue(null, value, false);
      }
    }
  }

  private static class TopRecordDropper
  extends SimpleTagValueWrapper {
    private int depth = 0;

    public TopRecordDropper(TagValueListener delegate) {
      super(delegate);
    }

    public void startRecord()
    throws ParserException {
      if(depth > 0) {
        super.startRecord();
      }

      depth++;
    }

    public void endRecord()
    throws ParserException {
      depth--;

      if(depth > 0) {
        super.endRecord();
      }
    }
  }

  private static class SPFeatureTableListener
  extends SimpleTagValueWrapper {
    private Pattern pat = Pattern.compile("(\\w+)\\s+(\\d+)\\s+(\\d+)");
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
      System.out.println(depth + " " + tag + " " + val);
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

