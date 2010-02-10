// $Id: FontChooser.java,v 1.2 2009/10/26 23:29:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// The FontChooser class is in the Public Domain, the code may be used
// for any purpose. It is provided as is with no warranty.
//
// The FontChooser class is based on the JFontChooser class written
// by: James Bardsley (torasin@torasin.com)
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.Component;
import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

public class FontChooser extends JDialog implements ActionListener, ListSelectionListener {

    private static final String   BOLD_ITALIC       = "Bold Italic";
    private static final String   ITALIC            = "Italic";
    private static final String   BOLD              = "Bold";
    private static final String   REGULAR           = "Regular";
    private static final String   DEFAULT_FONT_NAME = "Sans";
    public static final long      serialVersionUID  = 62256323L;
    private static final String[] STYLE             = { REGULAR, BOLD, ITALIC, BOLD_ITALIC };
    private static final String[] SIZE              = { "3", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22",
            "24", "26", "28", "36", "72"           };
    private static final int      OK_OPTION         = 1;
    private static final int      CANCEL_OPTION     = 2;
    private Font                  _font;
    private int                   _option;
    private String                _type;
    private int                   _style;
    private int                   _size;
    private final JList           _font_list        = new JList( Util.getAvailableFontFamiliesSorted() );
    private final JList           _style_list       = new JList( STYLE );
    private final JList           _size_list        = new JList( SIZE );
    private final JTextField      _fonts_tf         = new JTextField();
    private final JTextField      _style_tf         = new JTextField();
    private final JTextField      _size_tf          = new JTextField();
    private final JLabel          _fonts_label      = new JLabel( "Font:" );
    private final JLabel          _style_label      = new JLabel( "Style:" );
    private final JLabel          _size_label       = new JLabel( "Size:" );
    private final JScrollPane     _font_jsp         = new JScrollPane( _font_list );
    private final JScrollPane     _style_jsp        = new JScrollPane( _style_list );
    private final JScrollPane     _size_jsp         = new JScrollPane( _size_list );
    private final JButton         _ok_button        = new JButton( "OK" );
    private final JButton         _cancel_button    = new JButton( "Cancel" );
    private final JTextField      _test_tf          = new JTextField( "AaBbZz012" );

    public FontChooser() {
        this( new Font( DEFAULT_FONT_NAME, Font.PLAIN, 12 ) );
    }

    public FontChooser( final Font font ) {
        final Container container = getContentPane();
        final JPanel panel = new JPanel();
        final TitledBorder panel_border = new TitledBorder( "Demo" );
        _font = font;
        _type = _font.getFontName();
        _style = _font.getStyle();
        _size = _font.getSize();
        _font_list.setSelectionMode( 0 );
        _style_list.setSelectionMode( 0 );
        _size_list.setSelectionMode( 0 );
        _font_jsp.setHorizontalScrollBarPolicy( ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        _style_jsp.setHorizontalScrollBarPolicy( ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        _size_jsp.setHorizontalScrollBarPolicy( ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER );
        panel.setBorder( panel_border );
        _fonts_tf.setBounds( 8, 32, 121, 20 );
        _font_jsp.setBounds( 8, 56, 121, 82 );
        _style_tf.setBounds( 136, 32, 121, 20 );
        _style_jsp.setBounds( 136, 56, 121, 82 );
        _size_tf.setBounds( 264, 32, 41, 20 );
        _size_jsp.setBounds( 264, 56, 41, 82 );
        _ok_button.setBounds( 320, 8, 89, 17 );
        _cancel_button.setBounds( 320, 32, 89, 17 );
        panel.setBounds( 320, 64, 89, 73 );
        container.add( _fonts_label );
        container.add( _fonts_tf );
        container.add( _font_jsp );
        container.add( _style_label );
        container.add( _style_tf );
        container.add( _style_jsp );
        container.add( _size_label );
        container.add( _size_tf );
        container.add( _size_jsp );
        container.add( _ok_button );
        container.add( _cancel_button );
        container.add( panel );
        _test_tf.setBounds( 8, 25, 73, 30 );
        panel.add( _test_tf );
        container.setLayout( null );
        panel.setLayout( null );
        setSize( 424, 177 );
        setResizable( false );
        setModal( true );
        _fonts_tf.addActionListener( this );
        _size_tf.addActionListener( this );
        _style_tf.addActionListener( this );
        _cancel_button.addActionListener( this );
        _ok_button.addActionListener( this );
        _font_list.addListSelectionListener( this );
        _style_list.addListSelectionListener( this );
        _size_list.addListSelectionListener( this );
    }

    public FontChooser( final String font_name, final int font_style, final int size ) {
        this( new Font( font_name, font_style, size ) );
    }

    public void actionPerformed( final ActionEvent e ) {
        if ( e.getSource() == _fonts_tf ) {
            boolean found = false;
            _type = _fonts_tf.getText();
            for( int i = 0; i < _font_list.getModel().getSize(); i++ ) {
                if ( ( ( String ) _font_list.getModel().getElementAt( i ) ).startsWith( _fonts_tf.getText().trim() ) ) {
                    _font_list.setSelectedIndex( i );
                    setScrollPos( _font_jsp, _font_list, i );
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                _font_list.clearSelection();
            }
            else {
                _test_tf.setFont( new Font( _type, _style, _size ) );
            }
        }
        else if ( e.getSource() == _size_tf ) {
            boolean found = false;
            parseSize();
            _test_tf.setFont( new Font( _type, _style, _size ) );
            for( int i = 0; i < _size_list.getModel().getSize(); i++ ) {
                if ( _size_tf.getText().trim().equals( _size_list.getModel().getElementAt( i ) ) ) {
                    _size_list.setSelectedIndex( i );
                    setScrollPos( _size_jsp, _size_list, i );
                    found = true;
                    break;
                }
            }
            if ( !found ) {
                _size_list.clearSelection();
            }
        }
        else if ( e.getSource() == _style_tf ) {
            if ( _style_tf.getText().equals( REGULAR ) ) {
                _style = Font.PLAIN;
            }
            else if ( _style_tf.getText().equals( BOLD ) ) {
                _style = Font.BOLD;
            }
            else if ( _style_tf.getText().equals( ITALIC ) ) {
                _style = Font.ITALIC;
            }
            else if ( _style_tf.getText().equals( BOLD_ITALIC ) ) {
                _style = Font.BOLD & Font.ITALIC;
            }
            _style_list.setSelectedIndex( _style );
            _test_tf.setFont( new Font( _type, _style, _size ) );
        }
        else if ( e.getSource() == _ok_button ) {
            parseSize();
            _option = OK_OPTION;
            _font = new Font( _type, _style, _size );
            setVisible( false );
        }
        else if ( e.getSource() == _cancel_button ) {
            _option = CANCEL_OPTION;
            setVisible( false );
        }
    }

    @Override
    public Font getFont() {
        return _font;
    }

    public String getFontName() {
        return _font.getFontName();
    }

    public int getFontSize() {
        return _font.getSize();
    }

    public int getFontStyle() {
        return _font.getStyle();
    }

    private void parseSize() {
        try {
            _size = ( Integer.parseInt( _size_tf.getText().trim() ) );
        }
        catch ( final Exception ex ) {
            // Ignore.
        }
        if ( _size < 1 ) {
            _size = 1;
        }
    }

    @Override
    public void setFont( final Font font ) {
        _font = font;
    }

    private void setScrollPos( final JScrollPane sp, final JList list, final int index ) {
        final int unit_size = sp.getVerticalScrollBar().getMaximum() / list.getModel().getSize();
        sp.getVerticalScrollBar().setValue( ( index - 2 ) * unit_size );
    }

    public int showDialog( final Component parent, final String title ) {
        boolean found = false;
        _option = CANCEL_OPTION;
        setTitle( title );
        _test_tf.setFont( new Font( _type, _style, _size ) );
        for( int i = 0; i < _font_list.getModel().getSize(); i++ ) {
            _font_list.setSelectedIndex( i );
            if ( _font.getFamily().equals( _font_list.getSelectedValue() ) ) {
                found = true;
                setScrollPos( _font_jsp, _font_list, i );
                break;
            }
        }
        if ( !found ) {
            _font_list.clearSelection();
        }
        _style_list.setSelectedIndex( _font.getStyle() );
        found = false;
        for( int i = 0; i < _size_list.getModel().getSize(); i++ ) {
            _size_list.setSelectedIndex( i );
            if ( _font.getSize() <= Integer.parseInt( ( String ) _size_list.getSelectedValue() ) ) {
                found = true;
                setScrollPos( _size_jsp, _size_list, i );
                break;
            }
        }
        if ( !found ) {
            _size_list.clearSelection();
        }
        setLocationRelativeTo( parent );
        setVisible( true );
        return _option;
    }

    public void valueChanged( final ListSelectionEvent e ) {
        if ( e.getSource() == _font_list ) {
            if ( _font_list.getSelectedValue() != null ) {
                _fonts_tf.setText( ( ( String ) ( _font_list.getSelectedValue() ) ) );
            }
            _type = _fonts_tf.getText();
            _test_tf.setFont( new Font( _type, _style, _size ) );
        }
        else if ( e.getSource() == _style_list ) {
            _style_tf.setText( ( ( String ) ( _style_list.getSelectedValue() ) ) );
            if ( _style_tf.getText().equals( REGULAR ) ) {
                _style = 0;
            }
            else if ( _style_tf.getText().equals( BOLD ) ) {
                _style = 1;
            }
            else if ( _style_tf.getText().equals( ITALIC ) ) {
                _style = 2;
            }
            else if ( _style_tf.getText().equals( BOLD_ITALIC ) ) {
                _style = 3;
            }
            _test_tf.setFont( new Font( _type, _style, _size ) );
        }
        else if ( e.getSource() == _size_list ) {
            if ( _size_list.getSelectedValue() != null ) {
                _size_tf.setText( ( ( String ) ( _size_list.getSelectedValue() ) ) );
            }
            _size = ( Integer.parseInt( _size_tf.getText().trim() ) );
            _test_tf.setFont( new Font( _type, _style, _size ) );
        }
    }
}
