package org.biojava.bio.program.formats;

import java.util.regex.Pattern;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.EcNumber;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.program.tagvalue.BoundaryFinder;
import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.MultiTagger;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.RegexSplitter;
import org.biojava.bio.program.tagvalue.SimpleTagValueWrapper;
import org.biojava.bio.program.tagvalue.TagDelegator;
import org.biojava.bio.program.tagvalue.TagValueContext;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ParserException;
import org.biojava.utils.lsid.LifeScienceIdentifier;

public class Ligand {
  private Ligand() {}

  public static final class Enzyme
  implements Format {
    private static final AnnotationType ANNO_TYPE;
    private static final LineSplitParser PARSER;
    private static final LifeScienceIdentifier LSID;

    static {
      LSID = LifeScienceIdentifier.valueOf(
        "open-bio.org", "format", "ligand/enzyme" );

      Location NONE = CardinalityConstraint.NONE;
      Location ANY = CardinalityConstraint.ANY;
      Location ONE = CardinalityConstraint.ONE;

      PropertyConstraint c_string =
        new PropertyConstraint.ByClass(String.class);
      PropertyConstraint c_ecNumber =
        new PropertyConstraint.ByClass(EcNumber.class);

      PARSER = new LineSplitParser(LineSplitParser.GENBANK);

      ANNO_TYPE = new AnnotationType.Impl();
      ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
      ANNO_TYPE.setConstraints("ENTRY", c_ecNumber, ONE);
      ANNO_TYPE.setConstraints("NAME", c_string, ONE);
      ANNO_TYPE.setConstraints("CLASS", c_string, ONE);
      ANNO_TYPE.setConstraints("SYSNAME", c_string, ANY);
      ANNO_TYPE.setConstraints("REACTION", c_string, ANY);
      ANNO_TYPE.setConstraints("SUBSTRATE", c_string, ANY);
      ANNO_TYPE.setConstraints("PRODUCT", c_string, ANY);
      ANNO_TYPE.setConstraints("COMMENT", c_string, ANY);
      ANNO_TYPE.setConstraints("REFERENCE", c_string, ANY);
      ANNO_TYPE.setConstraints("PATHWAY", c_string, ANY);
      ANNO_TYPE.setConstraints("GENES", c_string, ANY);
      ANNO_TYPE.setConstraints("DISEASE", c_string, ANY);
      ANNO_TYPE.setConstraints("MOTIF", c_string, ANY);
      ANNO_TYPE.setConstraints("STRUCTURES", c_string, ANY);
      ANNO_TYPE.setConstraints("DBLINKS", c_string, ANY);
    }

    public ParserListener getParserListener(TagValueListener listener) {
      ChangeTable changeTable = new ChangeTable();

      changeTable.setChanger("ENTRY", new ChangeTable.Changer() {
        public Object change(Object value) {
          String sv = (String) value;
          sv = sv.substring("EC ".length(), sv.length());
          int spc = sv.indexOf(" ");
          if(spc != -1) { // remove obsolete - should publish this better
            sv = sv.substring(0, spc);
          }
          return EcNumber.Impl.valueOf(sv);
        }
      });

      ValueChanger valueChanger = new ValueChanger(listener, changeTable);

      DollarStringCatter catter = new DollarStringCatter(valueChanger);

      final Pattern ref_start = Pattern.compile("\\d+(\\s+\\[[^\\]]\\])?");

      TagDelegator tagDelegator = new TagDelegator(valueChanger);
      tagDelegator.setListener("NAME", catter);
      tagDelegator.setListener("SYSNAME", catter);
      tagDelegator.setListener("REACTION", catter);
      tagDelegator.setListener("SUBSTRATE", catter);
      tagDelegator.setListener("PRODUCT", catter);
      tagDelegator.setListener("DEFINITION", catter);
      tagDelegator.setListener("REFERENCE", new MultiTagger(valueChanger,
        new BoundaryFinder() {
          public boolean dropBoundaryValues() { return true; }
          public boolean isBoundaryStart(Object value) {
            // we should realy extract the UI entries
            String sv = (String) value;
            return ref_start.matcher(sv).matches();
          }
          public boolean isBoundaryEnd(Object value) { return false; }
        })
      );
      return new ParserListener(PARSER, tagDelegator);
    }

    
    public AnnotationType getType() {
      return ANNO_TYPE;
    }

    public LifeScienceIdentifier getLSID() {
      return LSID;
    }
  }

  public static class Reaction
  implements Format {
    private static final AnnotationType ANNO_TYPE;
    private static final LineSplitParser PARSER;
    private static final LifeScienceIdentifier LSID;

    static {
      LSID = LifeScienceIdentifier.valueOf(
              "open-bio.org", "format", "ligand/reaction" );

      PARSER = new LineSplitParser(LineSplitParser.GENBANK);

      Location NONE = CardinalityConstraint.NONE;
      Location ANY = CardinalityConstraint.ANY;
      Location ONE = CardinalityConstraint.ONE;
      Location ONE_OR_MORE = CardinalityConstraint.ONE_OR_MORE;

      PropertyConstraint c_string =
        new PropertyConstraint.ByClass(String.class);
      PropertyConstraint c_ecNumber =
        new PropertyConstraint.ByClass(EcNumber.class);

      ANNO_TYPE = new AnnotationType.Impl();
      ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
      ANNO_TYPE.setConstraints("ENTRY", c_string, ONE);
      ANNO_TYPE.setConstraints("NAME", c_string, ANY);
      ANNO_TYPE.setConstraints("DEFINITION", c_string, ONE_OR_MORE);
      ANNO_TYPE.setConstraints("EQUATION", c_string, ONE_OR_MORE);
      ANNO_TYPE.setConstraints("PATHWAY", c_string, ANY);
      ANNO_TYPE.setConstraints("ENZYME", c_ecNumber, ANY);
    }

    public ParserListener getParserListener(TagValueListener listener) {
      ChangeTable changeTable = new ChangeTable();
      changeTable.setChanger("ENZYME", FormatTools.EC_FROM_STRING);

      ValueChanger valueChanger = new ValueChanger(listener, changeTable);

      ChangeTable ct2 = new ChangeTable();
      ct2.setSplitter("ENZYME", new RegexSplitter(Pattern.compile("\\S+"), 0));

      ValueChanger vc2 = new ValueChanger(valueChanger, ct2);

      return new ParserListener(PARSER, vc2);
    }

    public AnnotationType getType() {
      return ANNO_TYPE;
    }

    public LifeScienceIdentifier getLSID() {
      return LSID;
    }
  }

  public static class Compound
  implements Format {
    private static final AnnotationType ANNO_TYPE;
    private static final LineSplitParser PARSER;
    private static final LifeScienceIdentifier LSID;

    static {
      LSID = LifeScienceIdentifier.valueOf(
              "open-bio.org", "format", "ligand/compound" );

      PARSER = new LineSplitParser(LineSplitParser.GENBANK);

      Location NONE = CardinalityConstraint.NONE;
      Location ANY = CardinalityConstraint.ANY;
      Location ONE = CardinalityConstraint.ONE;
      
      PropertyConstraint c_string =
        new PropertyConstraint.ByClass(String.class);
      PropertyConstraint c_ecNumber =
        new PropertyConstraint.ByClass(EcNumber.class);

      ANNO_TYPE = new AnnotationType.Impl();
      ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
      ANNO_TYPE.setConstraints("ENTRY", c_string, ONE);
      ANNO_TYPE.setConstraints("NAME", c_string, ONE);
      ANNO_TYPE.setConstraints("FORMULA", c_string, ANY);
      ANNO_TYPE.setConstraints("PATHWAY", c_string, ANY);
      ANNO_TYPE.setConstraints("REACTION", c_string, ANY);
      ANNO_TYPE.setConstraints("ENZYME", c_ecNumber, ANY);
      ANNO_TYPE.setConstraints("STRUCTURES", c_string, ANY);
      ANNO_TYPE.setConstraints("DBLINKS", c_string, ANY);
    }

    public ParserListener getParserListener(TagValueListener listener) {
      ChangeTable changeTable = new ChangeTable();
      changeTable.setChanger("ENZYME", FormatTools.EC_FROM_STRING);

      ValueChanger valueChanger = new ValueChanger(listener, changeTable);

      ChangeTable ct2 = new ChangeTable();
      RegexSplitter wordSplitter = new RegexSplitter(Pattern.compile("\\S+"), 0);
      ct2.setSplitter("ENZYME", wordSplitter);
      ct2.setSplitter("REACTION", wordSplitter);

      ValueChanger vc2 = new ValueChanger(valueChanger, ct2);

      return new ParserListener(PARSER, vc2);
    }

    public AnnotationType getType() {
      return ANNO_TYPE;
    }

    public LifeScienceIdentifier getLSID() {
      return LSID;
    }
  }
  private static class DollarStringCatter
  extends SimpleTagValueWrapper {
    private StringBuffer sBuff;

    public DollarStringCatter(TagValueListener listener) {
      super(listener);
      sBuff = new StringBuffer();
    }

    public void startTag(Object tag)
    throws ParserException {
      super.startTag(tag);

      sBuff.setLength(0);
    }

    public void value(TagValueContext ctxt, java.lang.Object value)
    throws ParserException {
      String sv = (String) value;
      if(sv.startsWith("$")) {
        sBuff.append(sv.substring(1, sv.length()));
      } else {
        sBuff.append(' ');
        sBuff.append(sv);
      }
    }

    public void endTag()
    throws ParserException {
      super.value(null, sBuff.toString());
      super.endTag();

      sBuff.setLength(0);
    }
  }
}
