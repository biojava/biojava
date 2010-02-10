package org.biojava.bio.program.tagvalue;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.utils.ParserException;

/**
 * <p>
 * A ValueChanger.Changer that returns a specific match value using a regex
 * Pattern.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class RegexChanger
  implements
    ChangeTable.Changer
{
  private Pattern pattern;
  private int matchGroup;

  /**
   * Create a new RegexChanger with a pattern.
   *
   * @param pattern  the Pattern used to split values
   * @param matchGroup the group to pull out - use 0 to pull out the whole match
   */
  public RegexChanger(Pattern pattern, int matchGroup) {
    this.pattern = pattern;
    this.matchGroup = matchGroup;
  }

  public Object change(Object value)
  throws ParserException {
    try {
      Matcher matcher = pattern.matcher(value.toString());
      matcher.find();
      return matcher.group(matchGroup);
    } catch (IllegalStateException e) {
      throw new ParserException(
        "Could not match " + pattern.pattern() + " to " + value,  e
      );
    }
  }
}

