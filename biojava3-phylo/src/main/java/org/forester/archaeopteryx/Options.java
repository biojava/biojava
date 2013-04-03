// $Id: Options.java,v 1.28 2009/12/09 20:13:50 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2009 Christian M. Zmasek
// Copyright (C) 2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.Font;

import org.forester.util.ForesterUtil;

/*
 * This is to hold changeable options.
 */
final class Options {

    static final double             MIN_CONFIDENCE_DEFAULT = 0.0;
    private boolean                 _show_node_boxes;
    private boolean                 _show_branch_length_values;
    private boolean                 _internal_number_are_confidence_for_nh_parsing;
    private boolean                 _show_scale;
    private boolean                 _show_overview;
    private boolean                 _antialias_screen;
    private boolean                 _antialias_print;
    private boolean                 _graphics_export_visible_only;
    private int                     _print_size_x;
    private int                     _print_size_y;
    private double                  _min_confidence_value;
    private boolean                 _print_black_and_white;
    private boolean                 _print_using_actual_size;
    private boolean                 _graphics_export_using_actual_size;
    private PHYLOGENY_GRAPHICS_TYPE _phylogeny_graphics_type;
    private CLADOGRAM_TYPE          _cladogram_type;
    private OVERVIEW_PLACEMENT_TYPE _ov_placement;
    private NODE_LABEL_DIRECTION    _node_label_direction;
    private Font                    _base_font;
    private boolean                 _match_whole_terms_only;
    private boolean                 _search_case_sensitive;
    private float                   _print_line_width;
    private boolean                 _inverse_search_result;
    private double                  _scale_bar_length;
    private short                   _number_of_digits_after_comma_for_confidence_values;
    private short                   _number_of_digits_after_comma_for_branch_length_values;
    private boolean                 _nh_parsing_replace_underscores;
    private boolean                 _nh_parsing_extract_pfam_taxonomy_codes;
    private boolean                 _nh_parsing_ignore_quotes;
    private boolean                 _editable;
    private boolean                 _background_color_gradient;

    private Options() {
        init();
    }

    final Font getBaseFont() {
        return _base_font;
    }

    final CLADOGRAM_TYPE getCladogramType() {
        return _cladogram_type;
    }

    final double getMinConfidenceValue() {
        return _min_confidence_value;
    }

    final NODE_LABEL_DIRECTION getNodeLabelDirection() {
        return _node_label_direction;
    }

    final short getNumberOfDigitsAfterCommaForBranchLengthValues() {
        return _number_of_digits_after_comma_for_branch_length_values;
    }

    final short getNumberOfDigitsAfterCommaForConfidenceValues() {
        return _number_of_digits_after_comma_for_confidence_values;
    }

    final OVERVIEW_PLACEMENT_TYPE getOvPlacement() {
        return _ov_placement;
    }

    final PHYLOGENY_GRAPHICS_TYPE getPhylogenyGraphicsType() {
        return _phylogeny_graphics_type;
    }

    final float getPrintLineWidth() {
        return _print_line_width;
    }

    final int getPrintSizeX() {
        return _print_size_x;
    }

    final int getPrintSizeY() {
        return _print_size_y;
    }

    final double getScaleBarLength() {
        return _scale_bar_length;
    }

    final private void init() {
        _show_node_boxes = false;
        _show_branch_length_values = false;
        _internal_number_are_confidence_for_nh_parsing = false;
        _show_scale = false;
        _antialias_screen = true;
        _antialias_print = true;
        _graphics_export_visible_only = false;
        _editable = true;
        _background_color_gradient = false;
        if ( Util.isUsOrCanada() ) {
            _print_size_x = Constants.US_LETTER_SIZE_X;
            _print_size_y = Constants.US_LETTER_SIZE_Y;
        }
        else {
            _print_size_x = Constants.A4_SIZE_X;
            _print_size_y = Constants.A4_SIZE_Y;
        }
        _min_confidence_value = MIN_CONFIDENCE_DEFAULT;
        _print_black_and_white = false;
        _print_using_actual_size = false;
        _graphics_export_using_actual_size = true;
        _phylogeny_graphics_type = PHYLOGENY_GRAPHICS_TYPE.RECTANGULAR;
        _base_font = new Font( Configuration.getDefaultFontFamilyName(), Font.PLAIN, 10 );
        _match_whole_terms_only = false;
        _search_case_sensitive = false;
        _print_line_width = Constants.PDF_LINE_WIDTH_DEFAULT;
        _show_overview = true;
        _ov_placement = OVERVIEW_PLACEMENT_TYPE.UPPER_LEFT;
        _node_label_direction = NODE_LABEL_DIRECTION.HORIZONTAL;
        _inverse_search_result = false;
        _scale_bar_length = 0.0;
        _number_of_digits_after_comma_for_branch_length_values = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_BRANCH_LENGTH_VALUES_DEFAULT;
        _number_of_digits_after_comma_for_confidence_values = Constants.NUMBER_OF_DIGITS_AFTER_COMMA_FOR_CONFIDENCE_VALUES_DEFAULT;
        _nh_parsing_replace_underscores = false;
        _nh_parsing_extract_pfam_taxonomy_codes = false;
        _nh_parsing_ignore_quotes = Constants.NH_PARSING_IGNORE_QUOTES_DEFAULT;
        _cladogram_type = Constants.CLADOGRAM_TYPE_DEFAULT;
    }

    final boolean isAntialiasPrint() {
        return _antialias_print;
    }

    final boolean isAntialiasScreen() {
        return _antialias_screen;
    }

    final boolean isEditable() {
        return _editable;
    }

    final boolean isExtractPfamTaxonomyCodesInNhParsing() {
        return _nh_parsing_extract_pfam_taxonomy_codes;
    }

    final boolean isGraphicsExportUsingActualSize() {
        return _graphics_export_using_actual_size;
    }

    final boolean isGraphicsExportVisibleOnly() {
        return _graphics_export_visible_only;
    }

    final boolean isInternalNumberAreConfidenceForNhParsing() {
        return _internal_number_are_confidence_for_nh_parsing;
    }

    final boolean isInverseSearchResult() {
        return _inverse_search_result;
    }

    final boolean isMatchWholeTermsOnly() {
        return _match_whole_terms_only;
    }

    final boolean isNhParsingIgnoreQuotes() {
        return _nh_parsing_ignore_quotes;
    }

    final boolean isPrintBlackAndWhite() {
        return _print_black_and_white;
    }

    final boolean isPrintUsingActualSize() {
        return _print_using_actual_size;
    }

    final boolean isReplaceUnderscoresInNhParsing() {
        return _nh_parsing_replace_underscores;
    }

    final boolean isSearchCaseSensitive() {
        return _search_case_sensitive;
    }

    final boolean isShowBranchLengthValues() {
        return _show_branch_length_values;
    }

    final boolean isShowNodeBoxes() {
        return _show_node_boxes;
    }

    final boolean isShowOverview() {
        return _show_overview;
    }

    final boolean isShowScale() {
        return _show_scale;
    }

    final void setAntialiasPrint( final boolean antialias_print ) {
        _antialias_print = antialias_print;
    }

    final void setAntialiasScreen( final boolean antialias_screen ) {
        _antialias_screen = antialias_screen;
    }

    final void setBaseFont( final Font base_font ) {
        _base_font = base_font;
    }

    final void setCladogramType( final CLADOGRAM_TYPE cladogram_type ) {
        _cladogram_type = cladogram_type;
    }

    final void setEditable( final boolean editable ) {
        _editable = editable;
    }

    final void setExtractPfamTaxonomyCodesInNhParsing( final boolean nh_parsing_extract_pfam_taxonomy_codes ) {
        _nh_parsing_extract_pfam_taxonomy_codes = nh_parsing_extract_pfam_taxonomy_codes;
    }

    final void setGraphicsExportUsingActualSize( final boolean graphics_export_using_actual_size ) {
        _graphics_export_using_actual_size = graphics_export_using_actual_size;
        if ( !graphics_export_using_actual_size ) {
            setGraphicsExportVisibleOnly( false );
        }
    }

    final void setGraphicsExportVisibleOnly( final boolean graphics_export_visible_only ) {
        _graphics_export_visible_only = graphics_export_visible_only;
        if ( graphics_export_visible_only ) {
            setGraphicsExportUsingActualSize( true );
        }
    }

    final void setInternalNumberAreConfidenceForNhParsing( final boolean internal_number_are_confidence_for_nh_parsing ) {
        _internal_number_are_confidence_for_nh_parsing = internal_number_are_confidence_for_nh_parsing;
    }

    final void setInverseSearchResult( final boolean inverse_search_result ) {
        _inverse_search_result = inverse_search_result;
    }

    final void setMatchWholeTermsOnly( final boolean search_whole_words_only ) {
        _match_whole_terms_only = search_whole_words_only;
    }

    final void setMinConfidenceValue( final double min_confidence_value ) {
        _min_confidence_value = min_confidence_value;
    }

    final void setNhParsingIgnoreQuotes( final boolean nh_parsing_ignore_quotes ) {
        _nh_parsing_ignore_quotes = nh_parsing_ignore_quotes;
    }

    final void setNodeLabelDirection( final NODE_LABEL_DIRECTION node_label_direction ) {
        _node_label_direction = node_label_direction;
    }

    final private void setNumberOfDigitsAfterCommaForBranchLength( final short number_of_digits_after_comma_for_branch_length_values ) {
        _number_of_digits_after_comma_for_branch_length_values = number_of_digits_after_comma_for_branch_length_values;
    }

    final private void setNumberOfDigitsAfterCommaForConfidenceValues( final short number_of_digits_after_comma_for_confidence_values ) {
        _number_of_digits_after_comma_for_confidence_values = number_of_digits_after_comma_for_confidence_values;
    }

    final void setOvPlacement( final OVERVIEW_PLACEMENT_TYPE ov_placement ) {
        _ov_placement = ov_placement;
    }

    final void setPhylogenyGraphicsType( final PHYLOGENY_GRAPHICS_TYPE phylogeny_graphics_type ) {
        _phylogeny_graphics_type = phylogeny_graphics_type;
    }

    final void setPrintBlackAndWhite( final boolean print_black_and_white ) {
        _print_black_and_white = print_black_and_white;
    }

    final void setPrintLineWidth( final float print_line_width ) {
        _print_line_width = print_line_width;
    }

    final void setPrintSizeX( final int print_size_x ) {
        _print_size_x = print_size_x;
    }

    final void setPrintSizeY( final int print_size_y ) {
        _print_size_y = print_size_y;
    }

    final void setPrintUsingActualSize( final boolean print_using_actual_size ) {
        _print_using_actual_size = print_using_actual_size;
    }

    final void setReplaceUnderscoresInNhParsing( final boolean nh_parsing_replace_underscores ) {
        _nh_parsing_replace_underscores = nh_parsing_replace_underscores;
    }

    final void setScaleBarLength( final double scale_bar_length ) {
        _scale_bar_length = scale_bar_length;
    }

    final void setSearchCaseSensitive( final boolean search_case_sensitive ) {
        _search_case_sensitive = search_case_sensitive;
    }

    final void setShowBranchLengthValues( final boolean show_branch_length_values ) {
        _show_branch_length_values = show_branch_length_values;
    }

    final void setShowNodeBoxes( final boolean show_node_boxes ) {
        _show_node_boxes = show_node_boxes;
    }

    final void setShowOverview( final boolean show_overview ) {
        _show_overview = show_overview;
    }

    final void setShowScale( final boolean show_scale ) {
        _show_scale = show_scale;
    }

    final static Options createDefaultInstance() {
        return new Options();
    }

    final static Options createInstance( final Configuration configuration ) {
        final Options instance = createDefaultInstance();
        if ( configuration != null ) {
            instance.setAntialiasScreen( configuration.isAntialiasScreen() );
            instance.setShowScale( configuration.isShowScale() );
            instance.setShowBranchLengthValues( configuration.isShowBranchLengthValues() );
            instance.setShowOverview( configuration.isShowOverview() );
            instance.setCladogramType( configuration.getCladogramType() );
            instance.setOvPlacement( configuration.getOvPlacement() );
            instance.setPrintLineWidth( configuration.getPrintLineWidth() );
            instance.setNodeLabelDirection( configuration.getNodeLabelDirection() );
            instance.setBackgroundColorGradient( configuration.isBackgroundColorGradient() );
            if ( configuration.getNumberOfDigitsAfterCommaForBranchLengthValues() >= 0 ) {
                instance.setNumberOfDigitsAfterCommaForBranchLength( configuration
                        .getNumberOfDigitsAfterCommaForBranchLengthValues() );
            }
            if ( configuration.getNumberOfDigitsAfterCommaForConfidenceValues() >= 0 ) {
                instance.setNumberOfDigitsAfterCommaForConfidenceValues( configuration
                        .getNumberOfDigitsAfterCommaForConfidenceValues() );
            }
            instance.setExtractPfamTaxonomyCodesInNhParsing( configuration.isExtractPfamTaxonomyCodesInNhParsing() );
            instance.setReplaceUnderscoresInNhParsing( configuration.isReplaceUnderscoresInNhParsing() );
            instance.setInternalNumberAreConfidenceForNhParsing( configuration
                    .isInternalNumberAreConfidenceForNhParsing() );
            instance.setNhParsingIgnoreQuotes( configuration.isNhParsingIgnoreQuotes() );
            instance.setEditable( configuration.isEditable() );
            if ( configuration.getMinConfidenceValue() != MIN_CONFIDENCE_DEFAULT ) {
                instance.setMinConfidenceValue( configuration.getMinConfidenceValue() );
            }
            if ( configuration.getGraphicsExportX() > 0 ) {
                instance.setPrintSizeX( configuration.getGraphicsExportX() );
            }
            if ( configuration.getGraphicsExportY() > 0 ) {
                instance.setPrintSizeY( configuration.getGraphicsExportY() );
            }
            if ( configuration.getBaseFontSize() > 0 ) {
                instance.setBaseFont( instance.getBaseFont().deriveFont( ( float ) configuration.getBaseFontSize() ) );
            }
            if ( !ForesterUtil.isEmpty( configuration.getBaseFontFamilyName() ) ) {
                instance.setBaseFont( new Font( configuration.getBaseFontFamilyName(), Font.PLAIN, instance
                        .getBaseFont().getSize() ) );
            }
            if ( configuration.getPhylogenyGraphicsType() != null ) {
                instance.setPhylogenyGraphicsType( configuration.getPhylogenyGraphicsType() );
            }
        }
        return instance;
    }

    static enum CLADOGRAM_TYPE {
        NON_LINED_UP, EXT_NODE_SUM_DEP, TOTAL_NODE_SUM_DEP;
    }

    static enum NODE_LABEL_DIRECTION {
        HORIZONTAL, RADIAL;
    }

    static enum OVERVIEW_PLACEMENT_TYPE {
        UPPER_LEFT( "upper left" ),
        UPPER_RIGHT( "upper right" ),
        LOWER_LEFT( "lower left" ),
        LOWER_RIGHT( "lower right" );

        private final String _name;

        private OVERVIEW_PLACEMENT_TYPE( final String name ) {
            _name = name;
        }

        @Override
        public String toString() {
            return _name;
        }

        public String toTag() {
            return toString().replaceAll( " ", "_" );
        }
    }

    static enum PHYLOGENY_GRAPHICS_TYPE {
        RECTANGULAR, TRIANGULAR, EURO_STYLE, ROUNDED, CONVEX, CURVED, UNROOTED, CIRCULAR;
    }

    public void setBackgroundColorGradient( final boolean background_color_gradient ) {
        _background_color_gradient = background_color_gradient;
    }

    public boolean isBackgroundColorGradient() {
        return _background_color_gradient;
    }
}
