<?php
# Adapted from https://gist.github.com/alexkingorg/2158428
# Can be run like this : php -f sort-hex-colors-dark2light.php colors="#FFAC4E, #F9C804"
#
parse_str(implode('&', array_slice($argv, 1)), $_GET);

function cf_sort_hex_colors($colors) {
    $map = array(
        '0' => 0,
        '1' => 1,
        '2' => 2,
        '3' => 3,
        '4' => 4,
        '5' => 5,
        '6' => 6,
        '7' => 7,
        '8' => 8,
        '9' => 9,
        'a' => 10,
        'b' => 11,
        'c' => 12,
        'd' => 13,
        'e' => 14,
        'f' => 15,
    );
    $c = 0;
    $sorted = array();
    foreach ($colors as $color) {
        $color = strtolower(str_replace('#', '', $color));
        if (strlen($color) == 6) {
            $condensed = '';
            $i = 0;
            foreach (preg_split('//', $color, -1, PREG_SPLIT_NO_EMPTY) as $char) {
                if ($i % 2 == 0) {
                    $condensed .= $char;
                }
                $i++;
            }
            $color_str = $condensed;
        }
        $value = 0;
        foreach (preg_split('//', $color_str, -1, PREG_SPLIT_NO_EMPTY) as $char) {
            $value += intval($map[$char]);
        }
        $value = str_pad($value, 5, '0', STR_PAD_LEFT);
        $sorted['_'.$value.$c] = '#'.$color;
        $c++;
    }
    ksort($sorted);
#    print_r(array_values($sorted));
    print_r(strtoupper(implode(", ", $sorted)));
    echo PHP_EOL;
    return $sorted;
}
cf_sort_hex_colors(explode(", ",$_GET['colors']));

?>
