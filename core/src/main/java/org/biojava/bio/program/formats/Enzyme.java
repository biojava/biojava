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
import org.biojava.bio.program.tagvalue.RegexFieldFinder;
import org.biojava.bio.program.tagvalue.RegexSplitter;
import org.biojava.bio.program.tagvalue.TagDelegator;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.lsid.LifeScienceIdentifier;

public class Enzyme
implements Format {
  private static final AnnotationType ANNO_TYPE;
  private static final LineSplitParser PARSER;
  private static final LifeScienceIdentifier LSID;

  static {
    LSID = LifeScienceIdentifier.valueOf("open-bio.org", "format", "enzyme");

    Location NONE = CardinalityConstraint.NONE;
    Location ANY = CardinalityConstraint.ANY;
    Location ONE = CardinalityConstraint.ONE;

    PARSER = new LineSplitParser(LineSplitParser.EMBL);

    PropertyConstraint c_string = new PropertyConstraint.ByClass(String.class);
    PropertyConstraint c_ecNumber = new PropertyConstraint.ByClass(EcNumber.class);

    AnnotationType DI = new AnnotationType.Impl();
    DI.setDefaultConstraints(PropertyConstraint.NONE, NONE);
    DI.setConstraints("Disease_name", c_string, ONE);
    DI.setConstraints("MIM:Number", c_string, ONE);
    PropertyConstraint c_diType = new PropertyConstraint.ByAnnotationType(DI);

    ANNO_TYPE = new AnnotationType.Impl();
    ANNO_TYPE.setDefaultConstraints(PropertyConstraint.NONE, NONE);
    ANNO_TYPE.setConstraints("ID", c_ecNumber, ONE);
    ANNO_TYPE.setConstraints("DE", c_string, ONE);
    ANNO_TYPE.setConstraints("AN", c_string, ANY);
    ANNO_TYPE.setConstraints("CA", c_string, ANY);
    ANNO_TYPE.setConstraints("CF", c_string, ANY);
    ANNO_TYPE.setConstraints("CC", c_string, ANY);
    ANNO_TYPE.setConstraints("DI", c_diType, ANY);
    ANNO_TYPE.setConstraints("PR", c_string, ANY);
    ANNO_TYPE.setConstraints("DR", c_string, ANY);
  }

  public ParserListener getParserListener(TagValueListener listener) {
    ChangeTable.Changer trailingDotStripper = new ChangeTable.Changer() {
      public Object change(Object value) {
        String val = (String) value;
        if(val.endsWith(".")) {
          return val.substring(0, val.length() - 1);
        } else {
          return val;
        }
      }
    };

    ChangeTable changeTable = new ChangeTable();

    changeTable.setChanger("ID", new ChangeTable.Changer() {
      public Object change(Object value) {
        return EcNumber.Impl.valueOf((String) value);
      }
    });
    changeTable.setChanger("AN", trailingDotStripper);
    changeTable.setChanger("DE", trailingDotStripper);
    changeTable.setChanger("CA", trailingDotStripper);
    changeTable.setChanger("CF", trailingDotStripper);
    changeTable.setSplitter("DR", new RegexSplitter(
      Pattern.compile("\\S+,\\s*\\S+;"),
      0 ));

    ValueChanger valueChanger = new ValueChanger(listener, changeTable);

    MultiTagger dotMultiTag = new MultiTagger(
      valueChanger,
      new BoundaryFinder() {
        public boolean dropBoundaryValues() { return false; }
        public boolean isBoundaryStart(Object value) { return false; }
        public boolean isBoundaryEnd(Object value) {
          return ((String) value).endsWith(".");
        }
      }
    );

    MultiTagger commentMultiTag = new MultiTagger(
      valueChanger,
      new BoundaryFinder() {
        public boolean dropBoundaryValues() { return false; }
        public boolean isBoundaryStart(Object value) {
          return ((String) value).startsWith("-!-");
        }
        public boolean isBoundaryEnd(Object value) { return false; }
      }
    );

    TagDelegator tagDelegator = new TagDelegator(valueChanger);
    tagDelegator.setListener("AN", dotMultiTag);
    tagDelegator.setListener("CA", dotMultiTag);
    tagDelegator.setListener("CC", commentMultiTag);
    tagDelegator.setListener("DI", new RegexFieldFinder(
      valueChanger,
      Pattern.compile("([^;]+);\\s*MIM:\\s*(\\S+)\\."),
      new String[] { "Disease_name", "MIM:Number" },
      false ));


    return new ParserListener(PARSER, tagDelegator);
  }

  public AnnotationType getType() {
    return ANNO_TYPE;
  }

  public LifeScienceIdentifier getLSID() {
    return LSID;
  }
}
