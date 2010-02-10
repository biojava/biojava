package org.biojava.bio.program.tagvalue;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.utils.ParserException;

/**
 * <p>
 * A TagValueParser that splits a line based upon a regular expression. There
 * are configuration parameters analgous to those in LineSplitParser for
 * configuring parsing details.
 * </p>
 *
 * @author Matthew Pocock
 * @author Keith James (enabled empty line EOR)
 * @since 1.3
 */
public class RegexParser
  implements
    TagValueParser
{
  private Pattern pattern = null;
  
  private int tagGroup = -1;
  
  private int valueGroup = -1;
  
  private String endOfRecord = null;
  
  private boolean trimTag = false;
  
  private boolean trimValue = false;

  private boolean continueOnEmptyTag = false;
  
  private boolean mergeSameTag = false;
  
  private String tag;
  
  /**
   * Create a new RegexParser with all boolean values set to false.
   */
  public RegexParser() {}
  
  /** 
   * Set the Pattern used to split lines.
   *
   * @param pattern  the Pattern used to split lines
   */
  public void setPattern(Pattern pattern) {
    this.pattern = pattern;
  }
  
  /**
   * Get the Pattern currently used to split lines.
   *
   * @return the current Pattern
   */
  public Pattern getPattern() {
    return pattern;
  }
  
  /**
   * Set the group number that will match the tag.
   *
   * @param group the tag group number
   */
  public void setTagGroup(int group) {
    this.tagGroup = group;
  }
  
  /**
   * Get the group number that matches the tag.
   *
   * @return the tag group number
   */
  public int getTagGroup() {
    return tagGroup;
  }

  /**
   * Set the group number that will match the value.
   *
   * @param group the value group number
   */
  public void setValueGroup(int group) {
    this.valueGroup = group;
  }
  
  /**
   * Get the group number that matches the value.
   *
   * @return the value group number
   */
  public int getValueGroup() {
    return valueGroup;
  }
  
  /**
   * Set the explicit end-of-record string.
   *
   * @param endOfRecord  the new endOfRecord String
   */
  public void setEndOfRecord(String endOfRecord) {
    this.endOfRecord = endOfRecord;
  }
  
  /**
   * Get the explicit end-of-record string.
   *
   * @return  the current endOfRecord String
   */
  public String getEndOfRecord() {
    return endOfRecord;
  }
  
  /**
   * Enable trimming of the tag using String.trim().
   *
   * @param trimTag  true if tags should be trimmed, false otherwise
   */
  public void setTrimTag(boolean trimTag) {
    this.trimTag = trimTag;
  }
  
  /**
   * See if trimming of tags is enabled.
   *
   * @return true if tag trimming is enabled, false otherwise
   */
  public boolean getTrimTag() {
    return trimTag;
  }

  /**
   * Enable trimming of the value using String.trim().
   *
   * @param trimValue  true if values should be trimmed, false otherwise
   */
  public void setTrimValue(boolean trimValue) {
    this.trimValue = trimValue;
  }
  
  /**
   * See if trimming of values is enabled.
   *
   * @return true if value trimming is enabled, false otherwise
   */
  public boolean getTrimValue() {
    return trimValue;
  }

  /**
   * Decide whether to treat empty tags as continuations of the previous non
   * -empty tag.
   *
   * @param continueOnEmptyTag  true if empty tags should be replaced, false
   *        otherwise
   */
  public void setContinueOnEmptyTag(boolean continueOnEmptyTag) {
    this.continueOnEmptyTag = continueOnEmptyTag;
  }
  
  /**
   * Report whether empty tags will be treated as continuations of the last non
   * -empty tag.
   *
   * @return true if empty tags will be replaced, false otherwise
   */
  public boolean getContinueOnEmptyTag() {
    return continueOnEmptyTag;
  }
  
  /**
   * Decide if multiple examples of a single tag should be merged into a single
   * start/endTag pair with multiple values, or multiple start/endTag pairs each
   * with a single value.
   *
   * @param mergeSameTag  true if tags will be merged, false otherwise
   */
  public void setMergeSameTag(boolean mergeSameTag) {
    this.mergeSameTag = mergeSameTag;
  }
  
  /**
   * Report whether empty tags will be treated as continuations of the last non
   * -empty tag.
   *
   * @return true if tags will be merged, false otherwise
   */
  public boolean getMergeSameTag() {
    return mergeSameTag;
  }
  
  public TagValue parse(Object o)
  throws ParserException {
    String line = o.toString();
    
    // Use of the special value for the EOR marker allows a blank line
    // to be used to delimit records. Many file formats are like this.
    if (endOfRecord != null) {
        if (endOfRecord == TagValueParser.EMPTY_LINE_EOR) {
            if (line.equals(TagValueParser.EMPTY_LINE_EOR)) {
                return null;
            }
        }
        else
        {
            if (line.startsWith(endOfRecord)) {
                return null;
            }
        }
    }
    
    Matcher matcher = pattern.matcher(line);
    if(!matcher.find()) {
      throw new ParserException("Could not match " + pattern.pattern() + " to " + line);
    }
    String tag = matcher.group(tagGroup);
    if(trimTag) {
      tag = tag.trim();
    }
    
    String value = matcher.group(valueGroup);
    if(trimValue) {
      value = value.trim();
    }
    
    if(continueOnEmptyTag && (tag.length() == 0)) {
      return new TagValue(this.tag, value, false);
    } else if(mergeSameTag && tag.equals(this.tag)) {
      return new TagValue(tag, value, false);
    } else {
      return new TagValue(this.tag = tag, value, true);
    }
  }
}
