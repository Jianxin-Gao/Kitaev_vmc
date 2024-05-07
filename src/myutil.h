//
// Created by 高建鑫 on 2024/4/14.
//

#ifndef KITAEV_MYUTIL_H
#define KITAEV_MYUTIL_H

size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool ParserBondDimension(int argc, char *argv[],
                         std::vector<size_t> &D_set);
bool IsFileExist(const std::string &);

enum ExtensionDir {
    ROW,
    COL
};

#endif //KITAEV_MYUTIL_H
