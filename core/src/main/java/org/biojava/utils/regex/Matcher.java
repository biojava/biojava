

package org.biojava.utils.regex;

import org.biojava.bio.seq.io.SymbolListCharSequence;
import org.biojava.bio.symbol.SymbolList;

/**
 * This class is analogous to java.util.Matcher except that it works
 * on SymbolLists instead of Strings.  All coordinates are in the 1-based
 * coordinate system used by SymbolLists.
 *
 * @author David Huen
 * @since 1.4
 */
public class Matcher
{
    private org.biojava.utils.regex.Pattern pattern;
    private java.util.regex.Matcher matcher;
    private SymbolList sl;

    Matcher(org.biojava.utils.regex.Pattern pattern, SymbolList sl)
    {
        this.pattern = pattern;
        this.sl = sl;

        matcher = pattern.getPattern().matcher(new SymbolListCharSequence(sl));
    }

    /**
     * Returns the index of the last character matched, plus one.
     * @return The index of the last character matched, plus one.
     */
    public int end() { return matcher.end() + 1; }
    /**
     * Returns the index of the last Symbol, plus one, 
     * of the subsequence captured by the given group during the previous match operation.
     * <p>
     * Capturing groups are indexed from left to right, starting at one. 
     * Group zero denotes the entire pattern, so the expression m.end(0) is equivalent to m.end().
     * @param group The index of a capturing group in this matcher's pattern.
     * @return The index of the last Symbol captured by the group, plus one, 
     * or -1 if the match was successful but the group itself did not match anything.
     */
    public int end(int group) throws IndexOutOfBoundsException 
    {
        int pos = matcher.end(group);
        if (pos == -1)
            return pos;
        else
            return pos + 1;
    }

    /**
     * Attempts to find the next subsequence of the input sequence that matches the pattern.
     * <p>
     * This method starts at the beginning of the input sequence or, 
     * if a previous invocation of the method was successful and the matcher 
     * has not since been reset, at the first Symbol not matched by the previous match.
     * If the match succeeds then more information can be obtained via the start, end, and group methods. 
     * @return true if, and only if, a subsequence of the input sequence matches this matcher's pattern.
     */
    public boolean find() { return matcher.find(); }

    /**
     * Resets this matcher and then attempts to find the next subsequence 
     * of the input sequence that matches the pattern, starting at the specified index.
     * <p>
     * If the match succeeds then more information can be obtained via the start, 
     * end, and group methods, and subsequent invocations of the find() method 
     * will start at the first Symbol not matched by this match. 
     * @return true if, and only if, a subsequence of the input sequence 
     * starting at the given index matches this matcher's pattern.
     */
    public boolean find(int start) throws IndexOutOfBoundsException { return matcher.find(start - 1); }
    /**
     * Returns the input subsequence matched by the previous match.
     * <p>
     * For a matcher m with input sequence s, the expressions m.group() 
     * and s.substring(m.start(), m.end()) are equivalent.
     * Note that some patterns, for example a*, match the empty SymbolList. 
     * This method will return the empty string when the pattern successfully matches the empty string in the input. 
     * @return The (possibly empty) subsequence matched by the previous match, in SymbolList form.
     */
    public SymbolList group()
    {
        return sl.subList(start(), end() - 1);
    }

    /**
     * Returns the input subsequence captured by the given group during the previous match operation.
     * <p>
     * For a matcher m, input sequence s, and group index g, the expressions 
     * m.group(g) and s.substring(m.start(g), m.end(g)) are equivalent.
     * Capturing groups are indexed from left to right, starting at one. 
     * Group zero denotes the entire pattern, so the expression m.group(0) is equivalent to m.group().
     * If the match was successful but the group specified failed to match 
     * any part of the input sequence, then null is returned. 
     * Note that some groups, for example (a*), match the empty string. 
     * This method will return the empty string when such a group successfully matches the emtpy string in the input. 
     * @return The (possibly empty) subsequence captured by the group during the previous match, 
     * or null if the group failed to match part of the input.
     */
    public SymbolList group(int group)
        throws IndexOutOfBoundsException
    {
        int start = matcher.start(group);
        int end = matcher.end(group);
        if ((start == -1) && (end == -1)) return null;
        else
            return sl.subList(start(group), end(group) - 1);
    }

    /**
     * Returns the number of capturing groups in this matcher's pattern.
     * <p>
     * Any non-negative integer smaller than the value returned 
     * by this method is guaranteed to be a valid group index for this matcher. 
     * @return The number of capturing groups in this matcher's pattern.
     */
    public int groupCount() { return matcher.groupCount(); }

    /**
     * Attempts to match the input SymbolList, starting at the beginning, against the pattern.
     * <p>
     * Like the matches method, this method always starts at the 
     * beginning of the input sequence; unlike that method, 
     * it does not require that the entire input sequence be matched.
     * If the match succeeds then more information can be obtained via the start, end, and group methods.
     * @return true if, and only if, a prefix of the input sequence matches this matcher's pattern.
     */
    public boolean lookingAt() { return matcher.lookingAt(); }

    /**
     * Attempts to match the entire input sequence against the pattern.
     * <p>
     * If the match succeeds then more information can be obtained via the start, end, and group methods. 
     * @return true if, and only if, the entire input sequence matches this matcher's pattern.
     */
    public boolean matches() { return matcher.matches(); }

    /**
     * Returns the Pattern object that compiled this Matcher.
     */
    public org.biojava.utils.regex.Pattern pattern()
    {
        return pattern;
    }

    /**
     * Resets this matcher.
     * <p>
     * Resetting a matcher discards all of its explicit state information and sets its append position to zero. 
     * @return this matcher.
     */
    public org.biojava.utils.regex.Matcher reset()
    {
        matcher = matcher.reset();
        return this;
    }

    /**
     * Resets this matcher with a new input SymbolList.
     * <p>
     * Resetting a matcher discards all of its explicit state information and sets its append position to zero. 
     * @return this matcher.
     */
    public org.biojava.utils.regex.Matcher reset(SymbolList sl)
    {
        this.sl = sl;
        matcher = matcher.reset(new SymbolListCharSequence(sl));
        return this;
    }

    /**
     * Returns the start index of the previous match.
     * @return The index of the first Symbol matched.
     */
    public int start() { return matcher.start() + 1; }
    /**
     * Returns the start index of the subsequence captured by the given group during the previous match operation.
     * <p>
     * Capturing groups are indexed from left to right, starting at one. 
     * Group zero denotes the entire pattern, so the expression m.start(0) is equivalent to m.start(). 
     * @param group The index of a capturing group in this matcher's pattern.
     * @return The index of the first character captured by the group, or -1 if the match was successful 
     * but the group itself did not match anything.
     */
    public int start(int group) 
    {
        int pos = matcher.start(group);
        if (pos == -1)
            return pos;
        else
            return pos + 1;
    }

}

