package org.biojava.utils;

/**
 *
 *
 * @author Matthew Pocock
 */
public class RepeatedCharSequence
        implements CharSequence
{
  private int length;
  private char character;
  private StringBuffer sbuf;
  private String string;

  public RepeatedCharSequence(int length, char character)
  {
    this.length = length;
    this.character = character;
  }

  public RepeatedCharSequence()
  {
    this(0, ' ');
  }

  public int getLength()
  {
    return length;
  }

  public void setLength(int length)
  {
    if(length < 0) {
      throw new IllegalArgumentException("Length can't be negative: " + length);
    }

    // optimization
    if(sbuf != null) {
      if(length < this.length) {
        sbuf.setLength(length);
      } else {
        for(int i = this.length; i < length; i++) {
          sbuf.append(character);
        }
      }
    }

    if(this.length != length) {
      string = null;
    }

    this.length = length;
  }

  public char getCharacter()
  {
    return character;
  }

  public void setCharacter(char character)
  {
    this.character = character;
    flush();
  }

  public void flush()
  {
    sbuf = null;
    string = null;
  }

  public int length()
  {
    return length;
  }

  public char charAt(int index)
  {
    if(index < 0 || index >= length) {
      throw new IndexOutOfBoundsException(
              "Attempted to read from index " + index + " of " + length);
    }

    return character;
  }

  public CharSequence subSequence(int start, int end)
  {
    if(
            start < 0 ||
            start >= length ||
            end < start ||
            end > length)
    {
      throw new IndexOutOfBoundsException(
              "Illegal indexes: " + start + "," + end +
              " of sequence length " + length);
    }

    return new RepeatedCharSequence(end - start, character);
  }

  public String toString()
  {
    if(string == null) {
      if(sbuf == null) {
        sbuf = new StringBuffer(length);
        for(int i = 0; i < length; i++) {
          sbuf.append(character);
        }
      }

      string = sbuf.toString();
    }

    return string;
  }
}
