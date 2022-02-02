#ifndef _SHORT_STRING_H
#define _SHORT_STRING_H

static const int SHORT_STRING_LENGTH = 255;
static const int VERY_SHORT_STRING_LENGTH = 7;

typedef char     ShortString[SHORT_STRING_LENGTH + 1];
typedef char VeryShortString[VERY_SHORT_STRING_LENGTH + 1];

#endif // #ifndef _SHORT_STRING_H
